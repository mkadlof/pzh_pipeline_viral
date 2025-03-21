#!/usr/bin/env python3

import argparse
import csv
import json
import os
import subprocess
from warnings import warn


def get_nextclade_db_version(input: str):
    if not os.path.exists(input):
        warn(f"File {input} does not exist.")
        return "UNKNOWN"
    with open(input, 'r') as f:
        data = json.load(f)
    return data['meta']['updated']


def parse_nextclade_output_csv2json(input: str, input2: str, output: str, sequence_source: str):
    if not os.path.exists(input):
        warn(f"File {input} does not exist.")
        output_json = {"status": "nie",
                       "database_name": "Nextclade",
                       "error_message": f"File {input} does not exist."}
        with open(output, 'w') as f:
            json.dump(output_json, f)
    else:
        db_version = get_nextclade_db_version(input2)
        program_version = \
        subprocess.run(["nextclade", "--version"], capture_output=True, text=True).stdout.strip().split()[-1]

        with open(input, 'r') as f:
            reader = csv.DictReader(f, delimiter=';')
            data = list(reader)
        if data[0]['qc.overallStatus'] == 'good':
            qc_status = "pass"
        else:
            qc_status = "warning"
        output_json = {"status": "tak",
                       "database_name": "Nextclade",
                       "database_version": db_version,
                       "program_name": "Nextclade",
                       "program_version": program_version,
                       "variant_name": data[0]['clade'],
                       "variant_qc_status": qc_status,
                       "sequence_source": sequence_source}
        with open(output, 'w') as f:
            json.dump(output_json, f, indent=4)

    print(f"File {output} saved.")


def main():
    parser = argparse.ArgumentParser(description='Convert Nextstrain output CSV to JSON')
    parser.add_argument('input', type=str, help='Input CSV file')
    parser.add_argument('input2', type=str, help='Input file nextclade.auspice.json')
    parser.add_argument('output', type=str, help='Output JSON file')
    parser.add_argument('sequence_source', type=str, help='Segment id or string "full_genome"')
    args = parser.parse_args()

    parse_nextclade_output_csv2json(args.input, args.input2, args.output, args.sequence_source)


if __name__ == '__main__':
    main()
