Before sending a pull request, please ensure that:

- `make all pyfaunus usagetips unittests test` passes with no errors
- there are no compiler warnings
- no temporary files (.o, .pqr, state etc) are included in the commit
- naming policies for classes, functions, and variables adhere to:
  - `MyClass` (i.e. PascalCase)
  - `myFunction()` (i.e. camelCase)
  - `my_variable` (i.e. snake_case)
- source files use soft-indendation of four. To setup automatic
  code formatting of proposed changes using `clang-format`, install
  the provided git commit hook:

  ``` bash
  ./scripts/git-pre-commit-format install
  ```
