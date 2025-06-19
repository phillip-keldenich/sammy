import tarfile
from pathlib import Path


def extract_xml_instances_to(xml_instance_basedir, 
                             xml_target_dir,
                             name_stem):
    basepath = Path(xml_instance_basedir)
    targetpath = Path(xml_target_dir)
    for archive in basepath.glob("**/*.tar.*"):
        if "020_samples" in str(archive):
            continue
        with tarfile.open(archive, "r:*") as tf:
            try:
                datename = archive.name[:-7]
                expected_model = f"{datename}/model.xml"
                with tf.extractfile(expected_model) as filestream:
                    rawdata = filestream.read()
                with open(targetpath / f"{name_stem}-{datename}.xml", "wb") as out:
                    out.write(rawdata)
            except KeyError:
                print(f"WARNING: archive '{archive}' does not have member '{expected_model}'!")

