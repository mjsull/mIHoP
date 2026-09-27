# SPDX-FileCopyrightText: 2026-present QPHIRE-genomics <Mitchell.Sullivan@health.qld.gov.au>
#
# SPDX-License-Identifier: MIT

import re
from pathlib import Path

import typer
from snk_cli import CLI

mihop = CLI(Path(__file__).parent.parent)

CHUNK_SIZE = 1024 * 1024

@mihop.app.command()
def download(
    directory: str = typer.Option(".", "--directory", "-d"),
    threads: int = typer.Option(1, "--threads", "-t"),
):
    """
    Download mIHoP clinical meta database to local
    """
    import hashlib
    import json
    import urllib.error
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor

    from rich.console import Console, Group
    from rich.live import Live
    from rich.progress import (
        BarColumn,
        DownloadColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
        TransferSpeedColumn,
    )
    from rich.table import Column

    MANIFEST_URL = "https://mihop-db.s3.us-west-2.amazonaws.com/manifest.json"

    with urllib.request.urlopen(MANIFEST_URL) as response:
        manifest = json.load(response)

    files = manifest["files"]

    console = Console()
    directory_path = Path(directory)
    directory_path.mkdir(parents=True, exist_ok=True)

    label_column = Column(width=15, no_wrap=True)
    overall_progress = Progress(
        TextColumn("[bold]{task.description}", justify="right", table_column=label_column),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    )
    file_progress = Progress(
        TextColumn("{task.description}", justify="right", table_column=label_column),
        BarColumn(),
        DownloadColumn(),
        TransferSpeedColumn(),
        console=console,
    )
    overall_task = overall_progress.add_task("Overall", total=len(files))

    def download_one(file: dict) -> None:
        filename = file["path"].split("/")[-1]
        expected_size = file["size_bytes"]
        dest = directory_path / filename
        tmp = dest.with_name(dest.name + ".part")

        if dest.exists() and dest.stat().st_size == expected_size:
            console.print(f"[green]{filename} already exists, skipping[/green]")
            overall_progress.advance(overall_task)
            return

        resume_from = 0
        hasher = hashlib.sha256()
        if tmp.exists():
            resume_from = tmp.stat().st_size
            if resume_from > expected_size:
                tmp.unlink()
                resume_from = 0
            else:
                with tmp.open("rb") as f:
                    hasher = hashlib.file_digest(f, "sha256")

        task = file_progress.add_task(filename, total=expected_size, completed=resume_from)

        if resume_from < expected_size:
            request = urllib.request.Request(file["url"])
            if resume_from:
                request.add_header("Range", f"bytes={resume_from}-")

            try:
                response = urllib.request.urlopen(request)
            except urllib.error.HTTPError as e:
                console.print(f"[red]{filename}: failed to download ({e}), skipping[/red]")
                file_progress.remove_task(task)
                overall_progress.advance(overall_task)
                return

            with response:
                if resume_from and response.status != 206:
                    resume_from = 0
                    hasher = hashlib.sha256()
                    file_progress.update(task, completed=0)

                with tmp.open("ab" if resume_from else "wb") as out:
                    while chunk := response.read(CHUNK_SIZE):
                        out.write(chunk)
                        hasher.update(chunk)
                        file_progress.advance(task, len(chunk))

        file_progress.remove_task(task)

        if hasher.hexdigest() != file["sha256"]:
            console.print(f"[red]{filename}: checksum mismatch, removing partial file[/red]")
            tmp.unlink(missing_ok=True)
        else:
            tmp.rename(dest)
            console.print(f"[green]Downloaded {filename}[/green]")
        overall_progress.advance(overall_task)

    with (
        Live(Group(overall_progress, file_progress), console=console),
        ThreadPoolExecutor(max_workers=max(threads, 1)) as pool,
    ):
        for future in [pool.submit(download_one, file) for file in files]:
            future.result()

ILLUMINA_READ_PATTERN = re.compile(
    r"^(?P<sample>.+)_S\d+(?:_L\d+)?_R(?P<read>[12])(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)

@mihop.app.command()
def prepare(
    fastq_dir: str, 
    output: str = typer.Option("samples.csv", "--output", "-o")
) -> Path:
    """
    Generate a samplesheet from FASTQ files found in `fastq_dir`.
    """
    fastq_dir_path, output_path = Path(fastq_dir), Path(output)
    if not fastq_dir_path.is_dir():
        raise FileNotFoundError(f"FASTQ directory does not exist: {fastq_dir_path}")

    suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    sep = "\t" if output_path.suffix.lower() == ".tsv" else ","

    samples: dict[str, dict[str, Path]] = {}
    for path in fastq_dir_path.rglob("*"):
        lower = path.name.lower()
        suffix = next((s for s in suffixes if lower.endswith(s)), None)
        if not path.is_file() or not suffix or lower.startswith("undetermined"):
            continue

        match = ILLUMINA_READ_PATTERN.match(path.name)
        sample, read = (match["sample"], match["read"]) if match else (path.name[: -len(suffix)], "1")

        reads = samples.setdefault(sample, {})
        if read in reads:
            raise ValueError(f"Multiple R{read} FASTQs found for sample {sample!r}: {reads[read]} and {path}")
        reads[read] = path

    with output_path.open("w", newline="") as out:
        out.write(sep.join(("sample", "fastq_1", "fastq_2")) + "\n")
        for sample, reads in sorted(samples.items()):
            out.write(sep.join((sample, str(reads.get("1", "")), str(reads.get("2", "")))) + "\n")

    return output_path