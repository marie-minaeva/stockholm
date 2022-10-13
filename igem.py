from flask import Flask, request
from run_container import run
from collections import defaultdict
from flask_cors import CORS, cross_origin

app = Flask(__name__)
app.config['CORS_HEADERS'] = 'Content-Type'

CORS(app)


def get_args(request):
    """
    Parsing input arguments to human-friendly format
    Parameters
    ----------
    request : object
        JSON request received from web-form
    Returns
    -------
    output : defaultdict
        Parsed arguments
    """
    output = defaultdict(str)
    for item in request.values.items():
        output[item[0]] = item[1]
    return output


# Worker app
@app.route("/igem-software", methods=['POST'])
@cross_origin(origin='*', headers=['Content-Type','Authorization'])
def igem_form():
    """
    Main flask app
    Returns
    -------
    data : defaultdict
        Results in JSON format
    """
    # Retrieves data
    args = get_args(request)
    type_inp = request.values["type_inp"]
    fasta_sequence = request.values["fasta_sequence"]
    pos = request.values["pos"]
    if pos == "":
        pos = "0"
    mandatory_mutation = request.values["mutation_mandatory"]
    number_of_mutant = request.values["number_of_mutant"]
    matrix = request.values["matrix"]
    preserve = request.values["preserve"]
    database = request.values["database"]
    # Call calculation tool
    data = run(type_inp=type_inp, fasta_sequence=fasta_sequence, pos=pos, number_of_mutant=number_of_mutant,
                   matrix=matrix, preserve=preserve, database=database,  mandatory_mutation=mandatory_mutation)

    data["args"] = args
    return data