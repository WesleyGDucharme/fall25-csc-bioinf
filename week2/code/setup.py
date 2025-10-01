from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os, shutil, sys, sysconfig

class CodonExtension(Extension):
    def __init__(self, name: str, source: str):
        self.source = source
        super().__init__(name, sources=[], language="c")

class BuildCodonExt(build_ext):
    def run(self):
        super().run()
        for ext in self.extensions:
            self.build_codon(ext)

    def build_codon(self, ext: CodonExtension):
        codon = shutil.which("codon")
        if not codon:
            raise RuntimeError("codon not found on PATH")

        out_so = Path(self.get_ext_fullpath(ext.name))
        out_so.parent.mkdir(parents=True, exist_ok=True)
        obj_path = out_so.with_suffix(out_so.suffix + ".o")

        # IMPORTANT: compile with PIC
        self.spawn([
            codon, "build", "--release", "--pyext", "--relocation-model=pic",
            "-module", ext.name,
            "-o", str(obj_path),
            ext.source
        ])

        codon_lib = Path.home() / ".codon" / "lib" / "codon"
        libraries = ["codonrt"]
        library_dirs = [str(codon_lib)]
        runtime_library_dirs = [str(codon_lib)]

        # If your linker needs explicit libpython, uncomment these:
        pyver = "python" + sysconfig.get_python_version()  # e.g., python3.12
        libraries.append(pyver)
        library_dirs.append((Path(sys.executable).resolve().parent.parent / "lib").as_posix())

        self.compiler.link_shared_object(
            [str(obj_path)],
            str(out_so),
            libraries=libraries,
            library_dirs=library_dirs,
            runtime_library_dirs=runtime_library_dirs if os.name != "nt" else None,
            extra_preargs=["-fPIC"],  # be explicit at link step as well
        )

setup(
    name="trviz_codon",
    version="0.0.1",
    packages=["trviz_codon"],  # you have trviz_codon/__init__.py
    ext_modules=[
        CodonExtension("trviz_codon.utils", "trviz_codon/utils.codon"),
    ],
    cmdclass={"build_ext": BuildCodonExt},
)
