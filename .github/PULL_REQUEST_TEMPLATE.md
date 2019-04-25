Please check that:

- `make all pyfaunus usagetips unittests test` passes with no errors
- there are no compiler warnings
- no temporary files (.o, .pqr, state etc) are included in the commit
- naming policies for classes, functions, and variables adhere to:

  Style      | Elements
  ---------- | -------------------------
  PascalCase | classes, namespaces
  camelCase  | functions
  snake_case | variables

- source files use soft-indendation of four. To setup automatic
  code formatting of proposed changes using `clang-format`, install
  the provided git commit hook:

  ``` bash
  ./scripts/git-pre-commit-format install
  ```
