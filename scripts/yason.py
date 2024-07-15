#!/usr/bin/env python

import sys

if sys.version_info < (3, 0):
    sys.stderr.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import os
import io
import json
import sys
import argparse
import warnings

try:
    jsonschema = True
    from jsonschema import Draft201909Validator
    from jsonschema.exceptions import best_match
except ImportError:
    jsonschema = False

try:
    import jinja2
except ImportError:
    warnings.warn("warning: missing jinja2 module")
    jinja2 = None


try:
    import ruamel.yaml as yaml
except ImportError:
    print("ruame.yaml is required", file=sys.stderr)
    sys.exit(1)
else:

    def yaml_safe_load(stream):
        return yaml.YAML(typ="safe").load(stream)

    def yaml_safe_dump(data, stream=None, **kwargs):
        return yaml.YAML(typ="safe").dump(data, stream=stream, **kwargs)


try:
    pygments = True
    from pygments import highlight
    from pygments.lexers import JsonLexer, YamlLexer
    from pygments.formatters import NullFormatter, Terminal256Formatter
except ImportError:
    pygments = False

parser = argparse.ArgumentParser(description="Convert json to yaml and vice versa")
parser.add_argument(
    "--indent", "-i", dest="indent", default=4, type=int, help="text indentation"
)
parser.add_argument(
    "--json", "-j", dest="alwaysjson", action="store_true", help="always output to json"
)
parser.add_argument(
    "infile",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help="json/yaml input",
)

if pygments:
    parser.add_argument(
        "--color", "-c", dest="color", action="store_true", help="syntax colored output"
    )

args = parser.parse_args()

terminal_red = "\033[91m"
terminal_yellow = "\033[93m"
terminal_default = "\033[39m"

if pygments:
    if args.color:
        formatter = Terminal256Formatter
    else:
        formatter = NullFormatter


def human_readable_path(path):
    out = str()
    for i in path:
        if not isinstance(i, int):
            out += str(i) + " -> "
    return out


def eprint(*args, **kwargs):
    """print to stderr"""
    print(*args, file=sys.stderr, **kwargs)


def print_table(schema):
    """pretty print schema as markdown table"""
    properties = schema.get("properties", "")
    if isinstance(properties, dict):
        eprint(terminal_yellow + "Need help, my young apprentice?\n" + terminal_red)
        eprint("{:25} | {:7} | {:50}".format(25 * "-", 7 * "-", 50 * "-"))
        eprint("{:25} | {:7} | {:50}".format("Property", "Type", "Description"))
        eprint("{:25} | {:7} | {:50}".format(25 * "-", 7 * "-", 50 * "-"))
        for key, value in properties.items():
            required = key in schema.get("required", [""])
            _property = key + str("*" if required else "")
            if isinstance(value, dict):
                _description = value.get("description", "")
                if "type" in value:
                    _default = value.get("default", "")
                    if _default != "":
                        _property = _property + "=" + str(_default)
                    _type = value["type"]
                    if isinstance(_type, list):
                        _type = "/".join(_type)
                    eprint(
                        "{:25} | {:7} | {}".format(
                            "`" + _property + "`", _type, _description
                        )
                    )
                else:
                    eprint(
                        "{:25} | {:7} | {}".format(
                            "`" + _property + "`", "n/a", _description
                        )
                    )
        eprint(terminal_default)  # restore terminal color


def validate_input(instance):
    """JSON schema checker"""
    if jsonschema:
        pathname = os.path.dirname(sys.argv[0])  # location of yason.py
        for subdir in [
            "/../docs",
            "/../share/faunus",
        ]:  # schema file can be in src or installed
            schemafile = os.path.abspath(pathname + subdir) + "/schema.yml"
            if os.path.exists(schemafile):
                with open(schemafile, "r") as f:
                    _schema = yaml_safe_load(f)
                    error = best_match(
                        Draft201909Validator(_schema).iter_errors(instance)
                    )
                    if error is not None:
                        eprint(
                            "{}{}\n".format(
                                human_readable_path(error.path), error.message
                            )
                        )
                        print_table(error.schema)
                        sys.exit(1)
                break


try:  # ... to read json
    file_as_str = args.infile.read()
    if jinja2:
        # additional files can be used with {% include "file" %}
        dirs = [os.getcwd(), os.path.dirname(os.path.realpath(__file__)) + "/../top"]
        loader = jinja2.FileSystemLoader(dirs)
        env = jinja2.Environment(loader=loader)
        file_as_str = env.from_string(file_as_str).render()  # render jinja2

    file_as_dict = json.loads(file_as_str)
    if "mcloop" in file_as_dict or "version" in file_as_dict:
        validate_input(file_as_dict)
    if args.alwaysjson:
        if pygments:
            file_as_str = highlight(file_as_str, JsonLexer(), formatter())
        print(file_as_str)
    else:
        stream = io.StringIO()
        yaml_safe_dump(file_as_dict, stream=stream)
        if pygments:
            stream = highlight(stream.getvalue(), YamlLexer(), formatter())
            print(stream)
        else:
            print(stream.getvalue())
except json.decoder.JSONDecodeError:
    try:  # ... to read yaml
        file_as_dict = yaml_safe_load(file_as_str)
        if "mcloop" in file_as_dict or "version" in file_as_dict:
            validate_input(file_as_dict)
        stream = json.dumps(file_as_dict, indent=args.indent)
        if pygments:
            stream = highlight(stream, JsonLexer(), formatter())
        print(stream)
    except yaml.parser.ParserError as exception:
        eprint("input error: invalid json or yaml format")
        eprint(exception)
        sys.exit(1)
