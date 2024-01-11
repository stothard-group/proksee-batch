import json
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import datetime
from typing import List, Any, Dict
import random
import string

# Global variables (which should later be defined in a config file).
genetic_code: str = "11"
orientation: str = "+"

def genbank_to_cgview_json(genbank_file: str, json_file: str) -> None:
    genbank_records: List[SeqRecord] = list(SeqIO.parse(genbank_file, "genbank"))

    now: datetime.datetime = datetime.datetime.now()

    map_id: str = ''.join(random.choices(string.ascii_lowercase + string.digits, k=40))

    json_data: Dict[str, Any] = {
        "cgview": {
            "version": "1.0",
            "created": now.strftime("%Y-%m-%d %H:%M:%S"),
            "id": map_id,
            "name": os.path.basename(genbank_file.rsplit('.', 1)[0]),
            "geneticCode": genetic_code,
            "settings": {},
            "backbone": {},
            "ruler": {},
            "dividers": {},
            "annotation": {},
            "sequence": {
                "contigs": []
            },
            "captions": [],
            "legend": {
                "items": []
            },
            "features": [],
            "tracks": [
                {
                    "name": "Features",
                    "separateFeaturesBy": "strand",
                    "position": "both",
                    "dataType": "feature",
                    "dataMethod": "source",
                    "dataKeys": "genbank-features"
                }
            ]
        }
    }

    for record in genbank_records:
        contig_data: Dict[str, Any] = {
            "name": record.name,
            "orientation": orientation,
            "length": len(record.seq),
            "seq": str(record.seq),

        }
        json_data["cgview"]["sequence"]["contigs"].append(contig_data)

        for feature in record.features:
            if feature.type == "CDS":
                feature_data: Dict[str, Any] = {
                    "type": feature.type,
                    "name": feature.qualifiers["locus_tag"][0],
                    "start": int(feature.location.start) + 1,
                    "stop": int(feature.location.end),
                    "strand": feature.location.strand,
                    "source": "genbank-features",
                    "contig": record.name,
                    "legend": feature.type
                }
                json_data["cgview"]["features"].append(feature_data)

    with open(json_file, 'w') as outfile:
        json.dump(json_data, outfile, indent=4)
