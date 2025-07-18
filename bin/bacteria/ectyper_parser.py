
import sys
import click
import json


@click.command()
@click.option('-i', '--input_file', help='[INPUT] a path to an input file with ectyper results',
              type=click.Path(), required=True)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(status, input_file, output, error=""):
    if status != "tak":
        json_output = {"program_name": "ectyper",
                       "status": status,
                        "error_message": error}
    else:
        json_output = {"program_name": "ectyper",
                       "status": status}
        slownik = {}
        header, body = open(input_file).readlines()
        line = body.rstrip().split("\t")
        h_antygen = line[6]
        h_antigen_data = []
        for element in h_antygen.split("/"):
            h_antigen_data.append({"antigen_id": element})


        o_antigen = line[5]
        o_antigen_data = []
        for element in o_antigen.split("/"):
            o_antigen_data.append({"antigen_id" : element})
        serovar_antigen = line[7]
        json_output = {"program_name": "ectyper",
                       "status": status,
                       "serotype_name" : serovar_antigen,
                       "antigen_o_data" : o_antigen_data,
                       "antigen_h_data" : h_antigen_data
                       }

    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, ensure_ascii=False, indent = 4))

    return True

if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])

