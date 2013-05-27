#!/bin/bash
if test -r "$1"
then
  echo ":filetype plugin indent on
:set expandtab
:set shiftwidth=2
:set softtabstop=2
:0
=G
:x!" > .indent.vim
  vim -s .indent.vim $1
  rm .indent.vim
fi
