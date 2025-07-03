if __name__ != "__main__":
    raise ImportError("Not meant to be imported!")


from pathlib import Path
import lzma
import json

raw_outputs = Path(__file__).parent / "02_output"
instance_dir = Path(__file__).parent.parent.parent / "full_instances"
output_files = list(raw_outputs.glob("*.out.json.xz"))
out_data = []

for output_file in output_files:
    with lzma.open(output_file, "rt") as f:
        data = json.load(f)
    if "outcome" in data and data["outcome"] == 'error':
        raise RuntimeError(f"There was an error ({output_file})!")
    lb = data['lb']
    ub = data['ub']
    time = data['events'][-1]['time']
    mes = len(data['mutually_exclusive_set'])
    instance_name = data['instance_name']
    out_data.append({"lb": lb, "ub": ub, "time": time, "mes": mes, "instance": instance_name})

with lzma.open("04_summary_data.json.xz", "wt") as f:
    json.dump(out_data, f)

