#!/usr/bin/env python3

"""
SELECT50 example: MACE-OMOL calculations with ASE and geomeTRIC.

SETUP
  Use Python 3.11 or newer in a dedicated virtual environment:
    python -m venv .venv
  Activate it on Linux/macOS: source .venv/bin/activate
  Activate it in Windows PowerShell: .venv/Scripts/Activate.ps1
  Install PyTorch for your CPU/CUDA setup using https://pytorch.org/get-started/locally/
  Then install the calculation packages in that environment:
    python -m pip install ase geometric numpy mace-torch
  The installed MACE version must provide mace.calculators.mace_omol.
  Check it without loading a model:
    python -c "from mace.calculators import mace_omol"
  Model extra_large downloads the MACE-OMOL checkpoint on first use. For an
  existing checkpoint, pass --model /path/to/checkpoint.model instead.
  MACE documentation: https://mace-docs.readthedocs.io/en/latest/

RUN
  Copy this script and your XYZ coordinates (angstrom) as molecule.xyz into
  a fresh job directory. Run each job in its own directory: the helper module,
  trajectories, frequency cache and result files use fixed or input-based names.
  Set --charge and --mult to the molecule's total charge and multiplicity.
  The following command requests a tight TS optimization on a CUDA GPU:
    python opt_mace_omol.py --xyz molecule.xyz --calculator mace-omol --model extra_large --device cuda --opt-type ts --opt-convergence tight --charge 0 --mult 1 > optimization.log 2>&1
  For CPU execution, replace --device cuda with --device cpu.
  For minimum optimization, replace --opt-type ts with --opt-type min.
  The default TS Hessian setting is first+last; minimum optimization does not
  calculate frequencies unless a Hessian option or a separate frequency run is used.
  After optimization, calculate frequencies and imaginary-mode animations with:
    python opt_mace_omol.py --xyz molecule_optimized_structure.xyz --calculator mace-omol --model extra_large --device cuda --opt-type freq --charge 0 --mult 1 > frequencies.log 2>&1
  For a single-point energy/forces calculation, use --opt-type spe.
  Help: python opt_mace_omol.py --help
  Append --dry-run to inspect a command without a calculation. Optimization
  dry runs still require geometric-optimize on PATH and write the helper module.

OUTPUT AND REPRODUCIBILITY
  Successful optimization exports molecule_optimization_trajectory.xyz and
  molecule_optimized_structure.xyz. Frequency jobs write summaries, CSVs and
  mode animations. This driver does not calculate Gibbs thermal corrections.
  Save the input, exact command, logs, package versions (python -m pip freeze),
  and checkpoint version/hash alongside results. This is a supplied workflow
  example; it is not a record of the settings for every SELECT50 calculation.

OTHER BACKENDS AND DEFAULTS
  The supplied driver also supports UMA and AIMNet2. Always specify
  --calculator mace-omol for the examples above: the default is UMA.
  Optional UMA dependencies: python -m pip install fairchem-core huggingface_hub
  Optional AIMNet dependencies: python -m pip install "aimnet[ase,hf]"
  Follow the selected model's access/setup instructions for those backends.
  INI job files are supported through --job; CLI arguments take precedence.
  Torch compilation is disabled by default for compatibility; use
  --allow-torch-compile only with a supported environment.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Windows/Torch safety defaults
# ---------------------------------------------------------------------------
# AIMNet/nvalchemiops can trigger torch.compile/TorchInductor.
# On Windows this can fail with TritonMissing. Disable by default.
os.environ.setdefault("TORCH_COMPILE_DISABLE", "1")
os.environ.setdefault("TORCHDYNAMO_DISABLE", "1")


# ---------------------------------------------------------------------------
# Helper module written to the run directory
# ---------------------------------------------------------------------------
# geomeTRIC's ASE engine needs an importable Python calculator class.
# For SPE/FREQ, this script also imports the same helper directly.
HELPER_MODULE = r"""
from __future__ import annotations

import os

os.environ.setdefault("TORCH_COMPILE_DISABLE", "1")
os.environ.setdefault("TORCHDYNAMO_DISABLE", "1")

from ase.calculators.calculator import Calculator, all_changes


class AIMNetCalculator(Calculator):
    implemented_properties = ["energy", "forces", "free_energy", "charges", "stress"]

    def __init__(
        self,
        model_name="aimnet2-2025",
        charge=0,
        mult=1,
        prefer_cuda=True,
        **kwargs,
    ):
        super().__init__()

        try:
            import torch
            import torch._dynamo

            torch._dynamo.config.disable = True

            if prefer_cuda and torch.cuda.is_available():
                try:
                    torch.set_default_device("cuda")
                except Exception:
                    pass

                print(f"AIMNet: CUDA available: {torch.cuda.get_device_name(0)}")
            else:
                print("AIMNet: CUDA unavailable or not requested; using AIMNet/PyTorch default.")

        except Exception as exc:
            print(f"AIMNet: torch/CUDA check failed: {exc}")

        try:
            from aimnet.calculators import AIMNet2ASE
        except ImportError:
            try:
                from aimnet.calculators.aimnet2ase import AIMNet2ASE
            except ImportError as exc:
                raise ImportError(
                    "AIMNet is not installed in this Python environment.\n"
                    "Install with:\n"
                    "  python -m pip install \"aimnet[ase,hf]\"\n"
                ) from exc

        self.model_name = model_name
        self.charge = int(charge)
        self.mult = int(mult)

        # Some AIMNet2ASE versions do not accept device= directly.
        # CUDA is requested through torch.set_default_device("cuda") above.
        self.inner = AIMNet2ASE(
            model_name,
            charge=self.charge,
            mult=self.mult,
            **kwargs,
        )

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.inner.calculate(atoms, properties, system_changes)
        self.results = dict(self.inner.results)


class UMACalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        model_name="uma-s-1p2",
        task_name="omol",
        device="cuda",
        charge=0,
        mult=1,
        seed=1,
        **kwargs,
    ):
        super().__init__()

        try:
            import torch

            print("UMA: CUDA available:", torch.cuda.is_available())
            if torch.cuda.is_available():
                print("UMA: CUDA device:", torch.cuda.get_device_name(0))
        except Exception:
            pass

        try:
            from fairchem.core import pretrained_mlip, FAIRChemCalculator
        except ImportError as exc:
            raise ImportError(
                "FAIR-Chem/UMA is not installed in this Python environment.\n"
                "Install with:\n"
                "  python -m pip install fairchem-core huggingface_hub\n"
            ) from exc

        self.model_name = model_name
        self.task_name = task_name
        self.device = device
        self.charge = int(charge)
        self.mult = int(mult)
        self.seed = int(seed)

        # FAIR-Chem API compatibility:
        # Some versions accept seed= in get_predict_unit(), while some do not.
        # Try with seed first, then retry without it.
        try:
            self.predictor = pretrained_mlip.get_predict_unit(
                self.model_name,
                device=self.device,
                seed=self.seed,
            )
        except TypeError as exc:
            if "seed" in str(exc):
                print(
                    "UMA: this FAIR-Chem version does not accept seed= in "
                    "get_predict_unit(); retrying without seed."
                )
                self.predictor = pretrained_mlip.get_predict_unit(
                    self.model_name,
                    device=self.device,
                )
            else:
                raise

        self.inner = FAIRChemCalculator(
            self.predictor,
            task_name=self.task_name,
        )

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        # FAIR-Chem molecular tasks use atoms.info metadata.
        atoms.info["charge"] = self.charge
        atoms.info["spin"] = self.mult

        self.inner.calculate(atoms, properties, system_changes)
        self.results = dict(self.inner.results)


class MACEOMOLCalculator(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(
        self,
        model_name="extra_large",
        device="cuda",
        charge=0,
        mult=1,
        default_dtype="float64",
        **kwargs,
    ):
        super().__init__()

        try:
            import torch

            print("MACE-OMOL: CUDA available:", torch.cuda.is_available())
            if torch.cuda.is_available():
                print("MACE-OMOL: CUDA device:", torch.cuda.get_device_name(0))
        except Exception:
            pass

        try:
            from mace.calculators import mace_omol
        except ImportError as exc:
            raise ImportError(
                "MACE-OMOL is not installed in this Python environment.\n"
                "Install with:\n"
                "  python -m pip install mace-torch\n"
            ) from exc

        self.model_name = model_name
        self.device = device
        self.charge = float(charge)
        self.mult = float(mult)
        self.default_dtype = default_dtype

        # MACE-OMOL foundation models are loaded through the mace_omol helper.
        # The official examples use atoms.info["charge"] and atoms.info["spin"].
        self.inner = mace_omol(
            model=self.model_name,
            device=self.device,
            default_dtype=self.default_dtype,
            **kwargs,
        )

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        # MACE-OMOL uses total charge and spin/multiplicity as global inputs.
        atoms.info["charge"] = self.charge
        atoms.info["spin"] = self.mult

        self.inner.calculate(atoms, properties, system_changes)
        self.results = dict(self.inner.results)
"""


# ---------------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------------

def str_to_bool(value: str) -> bool:
    value = str(value).strip().lower()

    if value in ["1", "true", "yes", "y", "on"]:
        return True

    if value in ["0", "false", "no", "n", "off"]:
        return False

    raise ValueError(f"Cannot interpret boolean value: {value}")


def write_helper_module(run_dir: Path) -> Path:
    helper_path = run_dir / "select50_mlip_calcs.py"
    helper_path.write_text(HELPER_MODULE, encoding="utf-8")
    return helper_path


def load_helper_module(helper_path: Path):
    spec = importlib.util.spec_from_file_location("select50_mlip_calcs", helper_path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import helper module from {helper_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules["select50_mlip_calcs"] = module
    spec.loader.exec_module(module)
    return module


def find_geometric() -> str:
    geometric_exe = shutil.which("geometric-optimize")

    if geometric_exe is None:
        raise RuntimeError(
            "Could not find geometric-optimize in this environment.\n"
            "Install it with:\n"
            "  python -m pip install geometric\n"
        )

    return geometric_exe


def print_cuda_status() -> None:
    try:
        import torch

        print("PyTorch version:", torch.__version__)
        print("CUDA available:", torch.cuda.is_available())

        if torch.cuda.is_available():
            print("CUDA device:", torch.cuda.get_device_name(0))
        else:
            print("CUDA device: none detected")

    except Exception as exc:
        print(f"Could not check CUDA status because torch import failed: {exc}")


def command_to_string(cmd: list[str]) -> str:
    return " ".join(f'"{x}"' if " " in x else x for x in cmd)


def split_raw_converge(value: str | None) -> list[str]:
    """
    Allows advanced geomeTRIC convergence usage.

    Example:
      --converge "set GAU_TIGHT"

    becomes:
      --converge set GAU_TIGHT
    """

    if value is None:
        return []

    value = str(value).strip()

    if not value:
        return []

    return value.split()


def convergence_preset_to_geometric(value: str | None) -> list[str]:
    """
    Convert friendly convergence names into geomeTRIC --converge arguments.

    These map to geomeTRIC Gaussian-style convergence presets.

      loose      -> set GAU_LOOSE
      normal     -> set GAU
      tight      -> set GAU_TIGHT
      verytight  -> set GAU_VERYTIGHT

    The raw --converge option overrides this.
    """

    if value is None:
        return []

    key = str(value).strip().lower()

    mapping = {
        "loose": ["set", "GAU_LOOSE"],
        "normal": ["set", "GAU"],
        "default": ["set", "GAU"],
        "tight": ["set", "GAU_TIGHT"],
        "verytight": ["set", "GAU_VERYTIGHT"],
        "very_tight": ["set", "GAU_VERYTIGHT"],
        "very-tight": ["set", "GAU_VERYTIGHT"],
    }

    if key not in mapping:
        raise ValueError(
            f"Unknown --opt-convergence preset: {value}. "
            "Use loose, normal, tight, or verytight."
        )

    return mapping[key]


# ---------------------------------------------------------------------------
# Job-file support
# ---------------------------------------------------------------------------

def write_inline_constraints_from_job(
    config: configparser.ConfigParser,
    job_path: Path,
) -> Path | None:
    """
    Reads inline constraints from:

    [constraints]
    text =
        $freeze
        distance 1 2

    and writes them to a real geomeTRIC constraints file.
    """

    if not config.has_section("constraints"):
        return None

    if not config.has_option("constraints", "text"):
        return None

    constraints_text = config.get("constraints", "text").strip()

    if not constraints_text:
        return None

    constraints_path = job_path.parent / f"{job_path.stem}_constraints.txt"

    with constraints_path.open("w", encoding="utf-8") as f:
        f.write(constraints_text)
        f.write("\n")

    return constraints_path


def apply_job_file(args: argparse.Namespace) -> argparse.Namespace:
    """
    Apply an INI-style job file.

    Rule:
      command-line arguments override job-file arguments.

    To make that possible, most argparse defaults are None and final defaults
    are applied after the job file is read.
    """

    if args.job is None:
        return args

    job_path = Path(args.job).resolve()

    if not job_path.exists():
        raise FileNotFoundError(f"Job file not found: {job_path}")

    config = configparser.ConfigParser(interpolation=None)
    config.read(job_path)

    def has(section: str, key: str) -> bool:
        return config.has_section(section) and config.has_option(section, key)

    def get(section: str, key: str):
        return config.get(section, key)

    def set_if_missing(attr: str, section: str, key: str, cast):
        if getattr(args, attr) is not None:
            return

        if has(section, key):
            setattr(args, attr, cast(get(section, key)))

    # calculation section
    set_if_missing("calculator", "calculation", "calculator", str)
    set_if_missing("model", "calculation", "model", str)
    set_if_missing("opt_type", "calculation", "opt_type", str)
    set_if_missing("charge", "calculation", "charge", int)
    set_if_missing("mult", "calculation", "mult", int)
    set_if_missing("output_prefix", "calculation", "output_prefix", str)

    # Inline constraints are used unless an explicit CLI or file constraint exists.
    if args.constraints is None:
        inline_constraints_path = write_inline_constraints_from_job(config, job_path)
        if inline_constraints_path is not None:
            args.constraints = inline_constraints_path

    # External constraints file in the job file overrides inline constraints,
    # unless --constraints was explicitly supplied at the command line.
    if has("calculation", "constraints") and args.cli_constraints_was_set is False:
        constraints_value = get("calculation", "constraints").strip()

        if constraints_value:
            constraints_path = Path(constraints_value)

            if not constraints_path.is_absolute():
                constraints_path = job_path.parent / constraints_path

            args.constraints = constraints_path

    # UMA section
    set_if_missing("task", "uma", "task", str)
    set_if_missing("device", "uma", "device", str)
    set_if_missing("seed", "uma", "seed", int)

    # MACE-OMOL section
    set_if_missing("mace_dtype", "mace", "dtype", str)

    # AIMNet section
    if args.prefer_cuda is None and has("aimnet", "prefer_cuda"):
        args.prefer_cuda = str_to_bool(get("aimnet", "prefer_cuda"))

    # geomeTRIC section
    set_if_missing("coordsys", "geometric", "coordsys", str)
    set_if_missing("hessian", "geometric", "hessian", str)
    set_if_missing("maxiter", "geometric", "maxiter", int)
    set_if_missing("trust", "geometric", "trust", float)
    set_if_missing("tmax", "geometric", "tmax", float)
    set_if_missing("converge", "geometric", "converge", str)
    set_if_missing("opt_convergence", "geometric", "opt_convergence", str)

    # Frequency section
    set_if_missing("freq_delta", "freq", "freq_delta", float)
    set_if_missing("freq_nfree", "freq", "freq_nfree", int)
    set_if_missing("freq_name", "freq", "freq_name", str)
    set_if_missing("imag_threshold", "freq", "imag_threshold", float)
    set_if_missing("write_modes", "freq", "write_modes", str)
    set_if_missing("mode_amplitude", "freq", "mode_amplitude", float)
    set_if_missing("mode_images", "freq", "mode_images", int)

    if args.reuse_freq_cache is None and has("freq", "reuse_freq_cache"):
        args.reuse_freq_cache = str_to_bool(get("freq", "reuse_freq_cache"))

    # Output section
    set_if_missing("output_prefix", "output", "output_prefix", str)

    if args.write_optimized_xyz is None and has("output", "write_optimized_xyz"):
        args.write_optimized_xyz = str_to_bool(get("output", "write_optimized_xyz"))

    if args.write_optimization_trajectory is None and has("output", "write_optimization_trajectory"):
        args.write_optimization_trajectory = str_to_bool(get("output", "write_optimization_trajectory"))

    # Torch section
    if args.allow_torch_compile is None and has("torch", "allow_torch_compile"):
        args.allow_torch_compile = str_to_bool(get("torch", "allow_torch_compile"))

    args.job = job_path
    return args


def apply_final_defaults(args: argparse.Namespace) -> argparse.Namespace:
    if args.calculator is None:
        args.calculator = "uma"

    if args.opt_type is None:
        args.opt_type = "ts"

    if args.model is None:
        if args.calculator == "uma":
            args.model = "uma-s-1p2"
        elif args.calculator == "aimnet":
            args.model = "aimnet2-2025"
        elif args.calculator == "mace-omol":
            args.model = "extra_large"

    if args.task is None:
        args.task = "omol"

    if args.device is None:
        args.device = "cuda"

    if args.mace_dtype is None:
        args.mace_dtype = "float64"

    if args.prefer_cuda is None:
        args.prefer_cuda = True

    if args.charge is None:
        args.charge = 0

    if args.mult is None:
        args.mult = 1

    if args.coordsys is None:
        args.coordsys = "tric"

    if args.opt_convergence is None:
        args.opt_convergence = "normal"

    if args.hessian is None:
        args.hessian = "auto"

    if args.maxiter is None:
        args.maxiter = 1250

    if args.trust is None:
        args.trust = 0.01

    if args.tmax is None:
        args.tmax = 0.03

    if args.seed is None:
        args.seed = 1

    if args.freq_delta is None:
        args.freq_delta = 0.01

    if args.freq_nfree is None:
        args.freq_nfree = 2

    if args.imag_threshold is None:
        args.imag_threshold = 10.0

    if args.write_modes is None:
        args.write_modes = "imaginary"

    if args.mode_amplitude is None:
        args.mode_amplitude = 0.5

    if args.mode_images is None:
        args.mode_images = 21

    if args.reuse_freq_cache is None:
        args.reuse_freq_cache = False

    if args.write_optimized_xyz is None:
        args.write_optimized_xyz = True

    if args.write_optimization_trajectory is None:
        args.write_optimization_trajectory = True

    if args.allow_torch_compile is None:
        args.allow_torch_compile = False

    if args.dry_run is None:
        args.dry_run = False

    if args.hessian == "auto":
        if args.opt_type == "ts":
            args.hessian = "first+last"
        else:
            args.hessian = "none"

    return args


def validate_args(args: argparse.Namespace) -> None:
    if args.calculator not in ["uma", "aimnet", "mace-omol"]:
        raise ValueError("--calculator must be uma, aimnet, or mace-omol")

    if args.opt_type not in ["ts", "min", "spe", "freq"]:
        raise ValueError("--opt-type must be ts, min, spe, or freq")

    if args.device not in ["cuda", "cpu"]:
        raise ValueError("--device must be cuda or cpu")

    if args.mace_dtype not in ["float32", "float64"]:
        raise ValueError("--mace-dtype must be float32 or float64")

    if args.coordsys not in ["tric", "dlc", "hdlc", "cart"]:
        raise ValueError("--coordsys must be tric, dlc, hdlc, or cart")

    if args.hessian not in ["none", "first", "first+last"]:
        raise ValueError("--hessian must be auto, none, first, or first+last")

    if args.opt_type == "ts" and args.hessian == "none":
        raise ValueError("TS optimization requires --hessian first or --hessian first+last")

    if args.opt_convergence not in ["loose", "normal", "tight", "verytight"]:
        raise ValueError("--opt-convergence must be loose, normal, tight, or verytight")

    if args.freq_nfree not in [2, 4]:
        raise ValueError("--freq-nfree must be 2 or 4")

    if args.write_modes not in ["none", "imaginary", "all"]:
        raise ValueError("--write-modes must be none, imaginary, or all")

    if args.mode_images < 3:
        raise ValueError("--mode-images must be at least 3")


# ---------------------------------------------------------------------------
# Calculator creation
# ---------------------------------------------------------------------------

def build_ase_settings(args: argparse.Namespace) -> tuple[str, dict]:
    if args.calculator == "uma":
        ase_class = "select50_mlip_calcs.UMACalculator"
        ase_kwargs = {
            "model_name": args.model,
            "task_name": args.task,
            "device": args.device,
            "charge": args.charge,
            "mult": args.mult,
            "seed": args.seed,
        }
        return ase_class, ase_kwargs

    if args.calculator == "aimnet":
        ase_class = "select50_mlip_calcs.AIMNetCalculator"
        ase_kwargs = {
            "model_name": args.model,
            "charge": args.charge,
            "mult": args.mult,
            "prefer_cuda": args.prefer_cuda,
        }
        return ase_class, ase_kwargs

    if args.calculator == "mace-omol":
        ase_class = "select50_mlip_calcs.MACEOMOLCalculator"
        ase_kwargs = {
            "model_name": args.model,
            "device": args.device,
            "charge": args.charge,
            "mult": args.mult,
            "default_dtype": args.mace_dtype,
        }
        return ase_class, ase_kwargs

    raise ValueError(f"Unknown calculator: {args.calculator}")


def make_calculator(args: argparse.Namespace, helper_module):
    if args.calculator == "uma":
        return helper_module.UMACalculator(
            model_name=args.model,
            task_name=args.task,
            device=args.device,
            charge=args.charge,
            mult=args.mult,
            seed=args.seed,
        )

    if args.calculator == "aimnet":
        return helper_module.AIMNetCalculator(
            model_name=args.model,
            charge=args.charge,
            mult=args.mult,
            prefer_cuda=args.prefer_cuda,
        )

    if args.calculator == "mace-omol":
        return helper_module.MACEOMOLCalculator(
            model_name=args.model,
            device=args.device,
            charge=args.charge,
            mult=args.mult,
            default_dtype=args.mace_dtype,
        )

    raise ValueError(f"Unknown calculator: {args.calculator}")


# ---------------------------------------------------------------------------
# geomeTRIC optimization
# ---------------------------------------------------------------------------

def build_geometric_command(args: argparse.Namespace) -> list[str]:
    geometric_exe = find_geometric()
    ase_class, ase_kwargs = build_ase_settings(args)

    cmd = [
        geometric_exe,
        str(args.xyz),
        "--engine",
        "ase",
        "--ase-class",
        ase_class,
        "--ase-kwargs",
        json.dumps(ase_kwargs),
        "--coordsys",
        args.coordsys,
        "--maxiter",
        str(args.maxiter),
        "--trust",
        str(args.trust),
        "--tmax",
        str(args.tmax),
    ]

    # Convergence behavior:
    #   --converge is the raw advanced override.
    #   --opt-convergence is the simple loose/normal/tight/verytight preset.
    if args.converge is not None:
        raw_converge_args = split_raw_converge(args.converge)

        if raw_converge_args:
            cmd.extend(["--converge", *raw_converge_args])
    else:
        preset_converge_args = convergence_preset_to_geometric(args.opt_convergence)

        if preset_converge_args:
            cmd.extend(["--converge", *preset_converge_args])

    if args.opt_type == "ts":
        cmd.extend(
            [
                "--transition",
                "yes",
                "--hessian",
                args.hessian,
            ]
        )

    elif args.opt_type == "min":
        if args.hessian != "none":
            cmd.extend(["--hessian", args.hessian])

    else:
        raise ValueError(
            f"geomeTRIC command only applies to opt-type ts or min, not {args.opt_type}"
        )

    if args.constraints is not None:
        cmd.extend(["--constraints", str(args.constraints)])

    return cmd


def find_geometric_trajectory(args: argparse.Namespace, run_dir: Path) -> Path | None:
    """
    Find geomeTRIC's optimization trajectory XYZ.

    geomeTRIC usually writes:
      <input_stem>_optim.xyz

    Prefer the expected filename in the run directory, then the input directory.
    A fallback is accepted only when there is one unique *_optim.xyz file.
    """

    candidate_files: list[Path] = []

    expected_run_dir = run_dir / f"{args.xyz.stem}_optim.xyz"
    expected_input_dir = args.xyz.parent / f"{args.xyz.stem}_optim.xyz"

    if expected_run_dir.is_file():
        return expected_run_dir

    if expected_input_dir.is_file():
        return expected_input_dir

    search_dirs: list[Path] = []

    for directory in [run_dir, args.xyz.parent]:
        if directory.exists() and directory not in search_dirs:
            search_dirs.append(directory)

    for search_dir in search_dirs:
        candidate_files.extend(search_dir.glob("*_optim.xyz"))

    candidate_files = sorted({p.resolve() for p in candidate_files if p.is_file()})

    if not candidate_files:
        return None

    if len(candidate_files) > 1:
        raise ValueError("Ambiguous optimization trajectories: "
                         + ", ".join(str(p) for p in candidate_files))

    return candidate_files[0]


def write_geometric_xyz_outputs(args: argparse.Namespace, run_dir: Path) -> tuple[Path | None, Path | None]:
    """
    After a successful geomeTRIC optimization, preserve/copy the optimization
    trajectory and write the final optimized structure as a single-frame XYZ.

    Outputs:
      <prefix>_optimization_trajectory.xyz
      <prefix>_optimized_structure.xyz

    For input:
      molecule.xyz

    Default outputs:
      molecule_optimization_trajectory.xyz
      molecule_optimized_structure.xyz
    """

    if not args.write_optimized_xyz and not args.write_optimization_trajectory:
        print("")
        print("Skipping optimization XYZ output writing.")
        return None, None

    try:
        from ase.io import read, write
    except ImportError as exc:
        print("")
        print("Warning: ASE is unavailable, so optimization XYZ outputs could not be written.")
        print(f"Reason: {exc}")
        return None, None

    prefix = args.output_prefix or args.xyz.stem

    trajectory_out = run_dir / f"{prefix}_optimization_trajectory.xyz"
    final_out = run_dir / f"{prefix}_optimized_structure.xyz"

    optim_xyz = find_geometric_trajectory(args, run_dir)

    if optim_xyz is None:
        print("")
        print("Warning: could not find a geomeTRIC *_optim.xyz trajectory file.")
        print("No optimization trajectory or final optimized XYZ was written.")
        return None, None

    written_trajectory: Path | None = None
    written_final: Path | None = None

    try:
        if args.write_optimization_trajectory:
            if optim_xyz.resolve() != trajectory_out.resolve():
                shutil.copy2(optim_xyz, trajectory_out)
                written_trajectory = trajectory_out
            else:
                written_trajectory = optim_xyz

        if args.write_optimized_xyz:
            final_atoms = read(str(optim_xyz), index=-1)
            write(str(final_out), final_atoms, format="xyz")
            written_final = final_out

        print("")
        print("Wrote optimization XYZ outputs:")

        if written_trajectory is not None:
            print(f"  trajectory:       {written_trajectory}")

        if written_final is not None:
            print(f"  optimized single: {written_final}")

        print(f"  source trajectory:{optim_xyz}")

        return written_trajectory, written_final

    except Exception as exc:
        print("")
        print(f"Warning: found geomeTRIC trajectory {optim_xyz}, but could not process it.")
        print(f"Reason: {exc}")
        return written_trajectory, written_final


def run_geometric(args: argparse.Namespace, run_dir: Path, env: dict) -> int:
    cmd = build_geometric_command(args)

    print("geomeTRIC command:")
    print(command_to_string(cmd))
    print("")

    if args.dry_run:
        print("Dry run requested; exiting without running geomeTRIC.")
        return 0

    result = subprocess.run(cmd, cwd=run_dir, env=env)

    if result.returncode != 0:
        print("")
        print("Optimization failed or stopped before convergence.", file=sys.stderr)
        print("Check the geomeTRIC .log file in this directory.", file=sys.stderr)
        print("No final optimized single-structure XYZ was written because geomeTRIC did not finish successfully.")
        return result.returncode

    # Preserve the optimization trajectory and write a clean single-frame optimized structure.
    try:
        write_geometric_xyz_outputs(args, run_dir)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("")
    print("Done.")

    if args.opt_type == "ts":
        print(
            "For a valid TS, check the final Hessian/frequency section for exactly "
            "one imaginary frequency with the expected reaction-coordinate motion."
        )
    else:
        if args.hessian == "none":
            print("Minimum optimization completed. No frequency calculation was requested.")
        else:
            print(
                "Minimum optimization completed. Check the final Hessian/frequency "
                "section for zero imaginary frequencies."
            )

    return 0


# ---------------------------------------------------------------------------
# SPE
# ---------------------------------------------------------------------------

def run_single_point(args: argparse.Namespace, helper_module) -> int:
    import numpy as np
    from ase.io import read, write

    atoms = read(str(args.xyz))
    atoms.calc = make_calculator(args, helper_module)

    energy_ev = atoms.get_potential_energy()
    forces = atoms.get_forces()

    prefix = args.output_prefix or args.xyz.stem

    result_txt = Path(f"{prefix}_spe_results.txt")
    forces_xyz = Path(f"{prefix}_spe_forces.extxyz")
    forces_csv = Path(f"{prefix}_spe_forces.csv")

    max_force = float(np.max(np.linalg.norm(forces, axis=1)))
    rms_force = float(np.sqrt(np.mean(forces**2)))

    with result_txt.open("w", encoding="utf-8") as f:
        f.write("Single-point MLIP calculation\n")
        f.write("=============================\n")
        f.write(f"Input XYZ:     {args.xyz}\n")
        f.write(f"Calculator:    {args.calculator}\n")
        f.write(f"Model:         {args.model}\n")
        f.write(f"Charge:        {args.charge}\n")
        f.write(f"Multiplicity:  {args.mult}\n")

        if args.calculator == "uma":
            f.write(f"UMA task:      {args.task}\n")
            f.write(f"UMA device:    {args.device}\n")

        f.write("\n")
        f.write(f"Energy / eV:   {energy_ev:.12f}\n")
        f.write(f"Max force / eV A^-1: {max_force:.8f}\n")
        f.write(f"RMS force / eV A^-1: {rms_force:.8f}\n")

    atoms.arrays["forces"] = forces
    write(str(forces_xyz), atoms)

    with forces_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["atom_index", "element", "Fx_eV_per_A", "Fy_eV_per_A", "Fz_eV_per_A"])

        for i, atom in enumerate(atoms):
            fx, fy, fz = forces[i]
            writer.writerow([i + 1, atom.symbol, fx, fy, fz])

    print("")
    print("Single-point calculation complete.")
    print(f"  Energy / eV: {energy_ev:.12f}")
    print(f"  Max force / eV A^-1: {max_force:.8f}")
    print(f"  RMS force / eV A^-1: {rms_force:.8f}")
    print("")
    print("Wrote:")
    print(f"  {result_txt}")
    print(f"  {forces_xyz}")
    print(f"  {forces_csv}")

    return 0


# ---------------------------------------------------------------------------
# FREQ
# ---------------------------------------------------------------------------

def frequency_display_value(freq) -> float:
    """
    ASE may represent imaginary frequencies as complex values.
    Report 0+500j as -500 cm^-1 for chemistry-style display.
    """

    real = float(getattr(freq, "real", freq))
    imag = float(getattr(freq, "imag", 0.0))

    if abs(imag) > abs(real):
        return -abs(imag)

    return real


def write_mode_animation_xyz(
    atoms,
    mode,
    output_path: Path,
    amplitude: float = 0.5,
    nimages: int = 21,
) -> None:
    import numpy as np

    symbols = atoms.get_chemical_symbols()
    base_positions = atoms.get_positions()

    mode = np.asarray(mode, dtype=float)
    norms = np.linalg.norm(mode, axis=1)
    max_norm = float(np.max(norms))

    if max_norm <= 1.0e-12:
        raise ValueError("Mode displacement is essentially zero; cannot write animation.")

    mode = mode / max_norm

    with output_path.open("w", encoding="utf-8") as f:
        for frame in range(nimages):
            phase = math.sin(2.0 * math.pi * frame / max(nimages - 1, 1))
            positions = base_positions + amplitude * phase * mode

            f.write(f"{len(symbols)}\n")
            f.write(f"mode animation frame {frame + 1}/{nimages}\n")

            for symbol, xyz in zip(symbols, positions):
                x, y, z = xyz
                f.write(f"{symbol:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def run_frequency(args: argparse.Namespace, helper_module) -> int:
    from contextlib import redirect_stdout

    from ase.io import read
    from ase.vibrations import Vibrations

    atoms = read(str(args.xyz))
    atoms.calc = make_calculator(args, helper_module)

    prefix = args.output_prefix or args.xyz.stem
    freq_name = args.freq_name or f"{prefix}_mlip_freq"

    vib = Vibrations(
        atoms,
        name=freq_name,
        delta=args.freq_delta,
        nfree=args.freq_nfree,
    )

    if not args.reuse_freq_cache:
        print(f"Cleaning previous frequency cache for name: {freq_name}")
        try:
            vib.clean()
        except Exception:
            pass

    print("")
    print("Running finite-difference frequency calculation.")
    print(f"  Frequency name: {freq_name}")
    print(f"  Displacement:   {args.freq_delta} Angstrom")
    print(f"  nfree:          {args.freq_nfree}")
    print("")

    vib.run()

    summary_txt = Path(f"{freq_name}_summary.txt")
    freqs_csv = Path(f"{freq_name}_frequencies.csv")

    with summary_txt.open("w", encoding="utf-8") as f:
        try:
            vib.summary(log=f)
        except TypeError:
            with redirect_stdout(f):
                vib.summary()

    freqs = vib.get_frequencies()

    imaginary_modes: list[int] = []

    with freqs_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "mode_index_zero_based",
                "mode_number_one_based",
                "frequency_display_cm^-1",
                "raw_real_cm^-1",
                "raw_imag_cm^-1",
                "is_imaginary",
            ]
        )

        for i, freq in enumerate(freqs):
            real = float(getattr(freq, "real", freq))
            imag = float(getattr(freq, "imag", 0.0))
            display = frequency_display_value(freq)
            is_imag = display < -abs(args.imag_threshold)

            if is_imag:
                imaginary_modes.append(i)

            writer.writerow([i, i + 1, display, real, imag, is_imag])

    print("")
    print("Frequency calculation complete.")
    print(f"  Number of modes: {len(freqs)}")
    print(f"  Imaginary modes below -{abs(args.imag_threshold)} cm^-1: {len(imaginary_modes)}")

    if imaginary_modes:
        print("  Imaginary mode numbers:")
        print("   ", ", ".join(str(i + 1) for i in imaginary_modes))

    print("")
    print("Wrote:")
    print(f"  {summary_txt}")
    print(f"  {freqs_csv}")

    if args.write_modes == "none":
        modes_to_write = []
    elif args.write_modes == "imaginary":
        modes_to_write = imaginary_modes
    elif args.write_modes == "all":
        modes_to_write = list(range(len(freqs)))
    else:
        raise ValueError(f"Unknown --write-modes option: {args.write_modes}")

    if modes_to_write:
        mode_dir = Path(f"{freq_name}_mode_animations")
        mode_dir.mkdir(exist_ok=True)

        print("")
        print("Writing mode animation XYZ files:")

        for mode_index in modes_to_write:
            mode = vib.get_mode(mode_index)
            display = frequency_display_value(freqs[mode_index])
            mode_file = mode_dir / f"mode_{mode_index + 1:04d}_{display:.1f}cm-1.xyz"

            try:
                write_mode_animation_xyz(
                    atoms=atoms,
                    mode=mode,
                    output_path=mode_file,
                    amplitude=args.mode_amplitude,
                    nimages=args.mode_images,
                )
                print(f"  {mode_file}")
            except Exception as exc:
                print(f"  Could not write mode {mode_index + 1}: {exc}")

    print("")
    print("Note:")
    print("  This frequency calculation is on the MLIP potential-energy surface.")
    print("  For final mechanistic validation, run a DFT frequency calculation in ORCA.")

    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="opt_mace_omol.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Run MLIP calculations with geomeTRIC + ASE.

Calculation modes:
  --opt-type ts      Transition-state optimization through geomeTRIC.
  --opt-type min     Minimum optimization through geomeTRIC.
  --opt-type spe     Single-point energy/force calculation. No optimization.
  --opt-type freq    Finite-difference frequency calculation. No optimization.

Backends:
  --calculator uma        FAIR-Chem UMA. Good broad-coverage OMol option.
  --calculator aimnet     AIMNet2. Use only for supported elements; not Bi.
  --calculator mace-omol  MACE-OMOL through mace-torch. Broad molecular MLIP with charge/spin inputs.

Convergence presets:
  --opt-convergence loose      geomeTRIC GAU_LOOSE
  --opt-convergence normal     geomeTRIC GAU
  --opt-convergence tight      geomeTRIC GAU_TIGHT
  --opt-convergence verytight  geomeTRIC GAU_VERYTIGHT

Job files:
  Use --job mlip.job to read settings from an INI-style file.
  Command-line arguments override job-file arguments.

Defaults after job-file parsing:
  calculator       uma
  UMA model        uma-s-1p2
  AIMNet model     aimnet2-2025
  MACE-OMOL model  extra_large
  opt-type         ts
  device           cuda
  coordsys         tric
  opt-convergence  normal
  hessian          auto
                   -> first+last for TS
                   -> none for min, spe, freq

After successful ts/min optimization, this script writes:
  <prefix>_optimization_trajectory.xyz
  <prefix>_optimized_structure.xyz
""",
        epilog="""
Examples:

  Help:
    python opt_mace_omol.py -h

  Check CUDA:
    python opt_mace_omol.py --check-cuda-only

  UMA TS optimization, normal convergence:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type ts --charge 0 --mult 1

  UMA TS optimization, tight convergence:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type ts --opt-convergence tight --charge 0 --mult 1

  MACE-OMOL TS optimization, tight convergence:
    python opt_mace_omol.py --xyz molecule.xyz --calculator mace-omol --model extra_large --opt-type ts --opt-convergence tight --charge 0 --mult 1

  UMA TS optimization, loose convergence for quick cleanup:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type ts --opt-convergence loose --charge 0 --mult 1

  Debug TS run with initial Hessian only:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type ts --hessian first

  UMA minimum optimization:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type min --charge 0 --mult 1

  MACE-OMOL TS optimization:
    python opt_mace_omol.py --xyz molecule.xyz --calculator mace-omol --model extra_large --opt-type ts --charge 0 --mult 1

  MACE-OMOL single-point calculation:
    python opt_mace_omol.py --xyz molecule.xyz --calculator mace-omol --model extra_large --opt-type spe --charge 0 --mult 1

  UMA single-point calculation:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --model uma-s-1p2 --opt-type spe --charge 0 --mult 1

  UMA frequency calculation:
    python opt_mace_omol.py --xyz molecule_optimized_structure.xyz --calculator uma --model uma-s-1p2 --opt-type freq --charge 0 --mult 1

  Frequency calculation with all mode animations:
    python opt_mace_omol.py --xyz molecule_optimized_structure.xyz --calculator uma --model uma-s-1p2 --opt-type freq --write-modes all

  Use raw geomeTRIC convergence instead of the friendly preset:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --opt-type ts --converge "set GAU_TIGHT"

  Disable final optimized single-frame extraction:
    python opt_mace_omol.py --xyz molecule.xyz --calculator uma --opt-type ts --no-write-optimized-xyz
""",
    )

    parser.add_argument("--job", type=Path, default=None, help="Optional INI-style job file.")
    parser.add_argument("--xyz", type=Path, default=None, help="Input XYZ file.")

    parser.add_argument("--calculator", choices=["uma", "aimnet", "mace-omol"], default=None, help="MLIP backend.")
    parser.add_argument("--model", default=None, help="Model name.")
    parser.add_argument("--opt-type", choices=["ts", "min", "spe", "freq"], default=None, help="Calculation type.")
    parser.add_argument("--task", default=None, help="UMA task name. Usually omol.")
    parser.add_argument("--device", choices=["cuda", "cpu"], default=None, help="Device for UMA or MACE-OMOL.")
    parser.add_argument("--mace-dtype", choices=["float32", "float64"], default=None, help="MACE-OMOL default dtype. float64 is safer; float32 is faster.")
    parser.add_argument("--prefer-cuda", action=argparse.BooleanOptionalAction, default=None, help="AIMNet CUDA preference.")

    parser.add_argument("--charge", type=int, default=None, help="Total charge.")
    parser.add_argument("--mult", type=int, default=None, help="Spin multiplicity.")

    parser.add_argument("--coordsys", choices=["tric", "dlc", "hdlc", "cart"], default=None, help="geomeTRIC coordinate system.")
    parser.add_argument("--hessian", choices=["auto", "none", "first", "first+last"], default=None, help="Hessian behavior.")
    parser.add_argument("--maxiter", type=int, default=None, help="Maximum geomeTRIC steps.")
    parser.add_argument("--trust", type=float, default=None, help="Initial geomeTRIC trust radius.")
    parser.add_argument("--tmax", type=float, default=None, help="Maximum geomeTRIC trust radius.")

    parser.add_argument(
        "--opt-convergence",
        choices=["loose", "normal", "tight", "verytight"],
        default=None,
        help="Friendly geomeTRIC optimization convergence preset.",
    )

    parser.add_argument(
        "--converge",
        default=None,
        help=(
            "Raw geomeTRIC convergence override. "
            "Example: --converge \"set GAU_TIGHT\". "
            "If set, this overrides --opt-convergence."
        ),
    )

    parser.add_argument("--constraints", type=Path, default=None, help="Optional external geomeTRIC constraints file.")

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="UMA seed. Automatically ignored if this FAIR-Chem version does not accept seed=.",
    )

    parser.add_argument("--freq-delta", type=float, default=None, help="Frequency displacement in Angstrom.")
    parser.add_argument("--freq-nfree", type=int, choices=[2, 4], default=None, help="Finite-difference stencil.")
    parser.add_argument("--freq-name", default=None, help="Base name for ASE frequency files.")
    parser.add_argument("--reuse-freq-cache", action=argparse.BooleanOptionalAction, default=None, help="Reuse ASE vibration cache.")
    parser.add_argument("--imag-threshold", type=float, default=None, help="Imaginary-mode threshold in cm^-1.")
    parser.add_argument("--write-modes", choices=["none", "imaginary", "all"], default=None, help="Write vibrational-mode animations.")
    parser.add_argument("--mode-amplitude", type=float, default=None, help="Mode animation amplitude in Angstrom.")
    parser.add_argument("--mode-images", type=int, default=None, help="Number of frames in mode animation.")

    parser.add_argument("--output-prefix", default=None, help="Output prefix for SPE/FREQ/optimization XYZ files.")

    parser.add_argument(
        "--write-optimized-xyz",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Write optimized single-frame XYZ after geomeTRIC optimization.",
    )

    parser.add_argument(
        "--write-optimization-trajectory",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Copy geomeTRIC optimization trajectory to <prefix>_optimization_trajectory.xyz.",
    )

    parser.add_argument("--allow-torch-compile", action=argparse.BooleanOptionalAction, default=None, help="Allow TorchDynamo/TorchInductor.")
    parser.add_argument("--check-cuda-only", action="store_true", help="Only print CUDA status and exit.")
    parser.add_argument("--dry-run", action=argparse.BooleanOptionalAction, default=None, help="Print command but do not run.")

    args = parser.parse_args()

    # Track whether user explicitly supplied --constraints.
    args.cli_constraints_was_set = args.constraints is not None

    return args


def main() -> int:
    args = parse_args()

    try:
        args = apply_job_file(args)
        args = apply_final_defaults(args)
        validate_args(args)
    except Exception as exc:
        print(f"ERROR in input/job configuration: {exc}", file=sys.stderr)
        return 1

    if args.allow_torch_compile:
        os.environ.pop("TORCH_COMPILE_DISABLE", None)
        os.environ.pop("TORCHDYNAMO_DISABLE", None)

    if args.check_cuda_only:
        print_cuda_status()
        return 0

    if args.xyz is None:
        print("ERROR: --xyz is required unless using --check-cuda-only.", file=sys.stderr)
        return 1

    args.xyz = args.xyz.resolve()

    if not args.xyz.exists():
        print(f"ERROR: XYZ file not found: {args.xyz}", file=sys.stderr)
        return 1

    if args.constraints is not None:
        args.constraints = Path(args.constraints).resolve()

        if not args.constraints.exists():
            print(f"ERROR: constraints file not found: {args.constraints}", file=sys.stderr)
            return 1

    if args.calculator == "aimnet":
        print(
            "\nWARNING: AIMNet2 standard models do not support Bi. "
            "Use --calculator uma or --calculator mace-omol for bismuth-containing systems.\n",
            file=sys.stderr,
        )

    print("")
    print("Requested calculation:")
    print(f"  input XYZ:   {args.xyz}")

    if args.job is not None:
        print(f"  job file:    {args.job}")

    print(f"  calculator:  {args.calculator}")
    print(f"  model:       {args.model}")
    print(f"  opt type:    {args.opt_type}")
    print(f"  charge:      {args.charge}")
    print(f"  mult:        {args.mult}")

    if args.opt_type in ["ts", "min"]:
        print(f"  hessian:     {args.hessian}")
        print(f"  coordsys:    {args.coordsys}")
        print(f"  convergence: {args.opt_convergence}")

        if args.converge is not None:
            print(f"  raw converge override: {args.converge}")

        print(f"  maxiter:     {args.maxiter}")
        print(f"  trust:       {args.trust}")
        print(f"  tmax:        {args.tmax}")
        print(f"  write trajectory XYZ: {args.write_optimization_trajectory}")
        print(f"  write final XYZ:      {args.write_optimized_xyz}")

        if args.constraints is not None:
            print(f"  constraints: {args.constraints}")

    if args.opt_type == "freq":
        print(f"  freq delta:  {args.freq_delta}")
        print(f"  freq nfree:  {args.freq_nfree}")
        print(f"  write modes: {args.write_modes}")

    if args.calculator == "uma":
        print(f"  UMA task:    {args.task}")
        print(f"  UMA device:  {args.device}")
    elif args.calculator == "mace-omol":
        print(f"  MACE device: {args.device}")
        print(f"  MACE dtype:  {args.mace_dtype}")
        print("  MACE-OMOL note: charge and spin are passed through atoms.info.")
    else:
        print(f"  AIMNet CUDA preference: {args.prefer_cuda}")
        print("  AIMNet note: this AIMNet2ASE version may not accept device= directly.")

    print("")
    print("CUDA status:")
    print_cuda_status()
    print("")

    run_dir = Path.cwd().resolve()
    helper_path = write_helper_module(run_dir)

    print("Wrote helper module:")
    print(f"  {helper_path}")
    print("")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(run_dir) + os.pathsep + env.get("PYTHONPATH", "")

    try:
        if args.opt_type in ["ts", "min"]:
            return run_geometric(args, run_dir, env)

        print("Direct ASE calculation selected; geomeTRIC will not be used.")
        print("")

        if args.dry_run:
            print("Dry run requested; exiting without running SPE/FREQ.")
            return 0

        helper_module = load_helper_module(helper_path)

        if args.opt_type == "spe":
            return run_single_point(args, helper_module)

        if args.opt_type == "freq":
            return run_frequency(args, helper_module)

        print(f"ERROR: unknown opt-type {args.opt_type}", file=sys.stderr)
        return 1

    except Exception as exc:
        print("")
        print(f"Calculation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
