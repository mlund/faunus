# Examples

To run the examples, the YAML input files must be piped through `yason.py`
to create the required JSON input. For examples:

~~~
yason.py water.yml | faunus
~~~

Similarly, JSON output can be converted to more readable YAML:

~~~
yason.py --color out.json
~~~

