# Examples

To run the examples, the YAML input files must be piped through `yason.py`
to create the required JSON input. For example,

    $ yason.py water.yml > water.json
    $ faunus -i water.json -o out.json

Or simply:

    $ yason.py water.yml | faunus

Reversely, JSON output can be converted to more readable YAML:

    $ yason.py out.json

