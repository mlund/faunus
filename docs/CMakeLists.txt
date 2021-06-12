# MANUAL

# Check for python packages required to build manual and parse yaml files
if (PYTHONINTERP_FOUND)
    # optional packages
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from bs4 import BeautifulSoup" RESULT_VARIABLE HAS_BS ERROR_QUIET)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pypandoc" RESULT_VARIABLE HAS_PYPANDOC ERROR_QUIET)

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from jsonschema import Draft7Validator" RESULT_VARIABLE HAS_JSONSCHEMA ERROR_QUIET)
    if (HAS_JSONSCHEMA EQUAL 1)
        message("python: missing `jsonschema` package is highly recommended")
    endif()

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ruamel_yaml as yaml" RESULT_VARIABLE HAS_YAML1 ERROR_QUIET)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ruamel.yaml as yaml" RESULT_VARIABLE HAS_YAML2 ERROR_QUIET)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import yaml" RESULT_VARIABLE HAS_YAML3 ERROR_QUIET)
    if ((HAS_YAML1 EQUAL 1) AND (HAS_YAML2 EQUAL 1) AND (HAS_YAML3 EQUAL 1))
        message("python: missing `ruamel_yaml` (preferred) or `yaml` packages are highly recommended")
    endif()
endif()

find_program(PANDOC pandoc)
if (PANDOC)
    if(NOT TEX_ENGINE)
        set(TEX_ENGINE xelatex CACHE STRING "TeX engine for PDF manual: xelatex lualatex")
    endif()
    set(TEX_ENGINE_OPTION --pdf-engine)
    execute_process(COMMAND ${PANDOC} ${TEX_ENGINE_OPTION}=${TEX_ENGINE}
        INPUT_FILE /dev/null ERROR_QUIET RESULT_VARIABLE PANDOC_CMD_SYNTAX)
    if(NOT PANDOC_CMD_SYNTAX EQUAL 0)
        # fallback to pandoc version 1 syntax
        set(TEX_ENGINE_OPTION --latex-engine)
    endif()
    add_custom_target(
        manual
        COMMENT "Building PDF manual"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs/_docs
        COMMAND cat
        header.md
        docs.md
        install.md
        running.md
        topology.md
        energy.md
        moves.md
        langevin.md
        analysis.md |
        grep -v ".notice--" |
        grep -v "Github-Improve_this_page" |
        sed 's.\\\\_._.g' | # this allows us to use \_ in markdown math
        ${PANDOC} --self-contained --highlight-style=pygments --template ../pandoc.tex
        -V date:${GIT_LATEST_TAG}""
        -N --toc ${TEX_ENGINE_OPTION}=${TEX_ENGINE} -o ${CMAKE_BINARY_DIR}/manual.pdf
        )
    set_target_properties(manual PROPERTIES EXCLUDE_FROM_ALL TRUE)

    add_custom_target(
        manual_html
        COMMENT "Building html manual"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs/_docs
        COMMAND cat
        header.md
        docs.md
        install.md
        running.md
        topology.md
        energy.md
        moves.md
        langevin.md
        analysis.md |
        grep -v ".notice--" |
        grep -v "Github-Improve_this_page" |
        sed 's.\\\\_._.g' | # this allows us to use \_ in markdown math
        ${PANDOC} --template ../GitHub.html5 -s --mathjax
        -M title:The\ Faunus\ User\ Guide
        -V date:${GIT_LATEST_TAG}""
        -N --toc -o ${CMAKE_BINARY_DIR}/manual.html
        )

    if (PYTHONINTERP_FOUND)
        if (HAS_BS EQUAL 0)
            if (HAS_PYPANDOC EQUAL 0)
                add_custom_target(
                    usagetips
                    ALL
                    COMMENT "extracting tips from manual"
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/extracttips.py)
                add_dependencies(usagetips manual_html)
            endif()
        endif()
    endif()
endif()

install(FILES
        ${CMAKE_BINARY_DIR}/manual.pdf
        ${CMAKE_BINARY_DIR}/manual.html
        ${CMAKE_SOURCE_DIR}/docs/schema.yml
        DESTINATION share/faunus OPTIONAL)
