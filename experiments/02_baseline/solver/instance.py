from pathlib import Path
from zipfile import ZipFile


def get_instance_from_archive(instance_name, out_dir,
                              archive_path):

    """
    Simple helper to parse instance
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with ZipFile(archive_path) as archive:
        # extract folder models/<instance_name> into out dir
        archive.extractall(
            path=out_dir,
            members=[
                f
                for f in archive.filelist
                if f.filename.startswith(f"models/{instance_name}")
            ],
        )

    out_path = Path(out_dir) / "models" / instance_name / "model.xml"

    if not out_path.exists():
        out_path = out_path.with_suffix(".dimacs")
        assert out_path.exists(), f"File {out_path} does not exist"

    return out_path
