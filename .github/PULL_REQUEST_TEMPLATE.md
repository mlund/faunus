Before sending a pull request, please ensure that:

- all tests pass, i.e. run `make test`
- there are no compiler warnings
- no temporary files (.o, .pqr, state etc) are included in the commit
- source files use soft-indendation of four. To setup automatic
  code formatting of proposed changes using `clang-format`, install
  the provided git commit hook:

  ``` bash
  ./scripts/git-pre-commit-format install
  ```
