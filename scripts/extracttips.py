#!/usr/bin/env python

#
# Converts all tables from manual.html to Markdown and stores
# them in a JSON file where the table text from the
# _upper left corner_ is used as keys.
# Long equations are filtered out.
#
# In addition, this adds a list of arbitrary ASCII to
# the `ascii` key.
#

import warnings
warnings.filterwarnings('ignore')
from bs4 import BeautifulSoup
import json
import pypandoc
import re

with open('manual.html', 'r') as f:
    data=f.read()
    soup = BeautifulSoup(data)

    headings = {}
    tables = {}  # database with markdown tables

    # http://www.ascii-art.de/ascii/
    # http://ascii.co.uk/art
    # http://www.oocities.org/spunk1111/wildlife.htm
    tables['ascii'] = [
    r'''
  ___________                               _________
   \_   _____/_____    __ __   ____   __ __ /   _____/
    |    __)  \__  \  |  |  \ /    \ |  |  \\_____  \
    |     \    / __ \_|  |  /|   |  \|  |  //        \
    \___  /   (____  /|____/ |___|  /|____//_______  /
        \/         \/             \/               \/''',

    r'''
    _                ___       _.--.
    \`.|\..----...-'`   `-._.-'_.-'`
    /  ' `         ,       __.--'
    )/' _/     \   `-_,   /
    `-'" `"\_  ,_.-;_.-\_ ',     fsc/as
        _.-'_./   {_.'   ; /
       {_.-``-'         {_/''',

    r"""
                            )/_
                  _.--..---"-,--c_
             \L..'           ._O__)_
     ,-.     _.+  _  \..--( /           a:f
       `\.-''__.-' \ (     \_
         `'''       `\__   /\
                     ')""",

    r'''
      ----------.................._____________  _  .-.
                                      _____.. . .   | |
                   _____....------""""             uuuuu
 ____....------""""                                |===|
                                                   |===|
                                                   |===|
                                                   |===|
 _a:f____________________________________________ .[__N]. _______''',

    r'''
               boing         boing         boing
     e-e           . - .         . - .         . - .
    (\_/)\       '       `.   ,'       `.   ,'       .
     `-'\ `--.___,         . .           . .          .
        '\( ,_.-'
           \\               "             "            a:f
           ^' ''',

    r'''
     "Mystery is the key to enchantment"             Armando Frazao
              ,d@@b,                    \.             a.k.a. Seal
     ._.__._._@@@@@@__...__.._..___._. _) `----._ ._..__._._.__.___._._
             -_-__-_-               _.'         e`.__
              -_-__-           _ ,-'..---~~~)/---'~~~
                _-       _ - '.,',',-       '
                          -  -_ -_ -''',

    r'''
                    ,.-----__
                 ,:::://///,:::-.
                /:''/////// ``:::`;/|/
               /'   ||||||     :://'`\
             .' ,   ||||||     `/(  e \
       -===~__-'\__X_`````\_____/~`-._ `.
                   ~~        ~~       `~-'  Seal''',

    r'''
                                ____
                               /\' .\    _____
                              /: \___\  / .  /\
   valkyrie                   \' / . / /____/..\
                               \/___/  \'  '\  /
                                        \'__'\/''',

    r'''
                                                  _  _
                                                 (\\( \
                                                  `.\-.)
                              _...._            _,-'   `-.
\                           ,'      `-._.---.,-'       .  \
 \`.                      ,'                               `.
  \ `-...__              /                           .   .:  y
   `._     ``--..__     /                           ,'`---._/
      `-._         ``--'                      |    /_
          `.._                   _            ;   <_ \
              `--.___             `.           `-._ \ \
                     `--<           `.     (\ _/)/ `.\/
                         \            \     `''',

    r'''
                  ,,__
        ..  ..   / o._)                   .---.
       /--'/--\  \-'||        .----.    .'     '.
      /        \_/ / |      .'      '..'         '-.
    .'\  \__\  __.'.'     .'          i-._
      )\ |  )\ |      _.'
     // \\ // \\
    ||_  \\|_  \\_
mrf '--' '--'' '--' ''',

    r'''
               __.....__
            .-'         '-.
          .'               '.
         /                   \
        /        |\           \
       ;        |V \_          ;
       |        |  ' \         ;
       ;        )   ,_\        |
       ;       /    |          ;
        \     /      \        /
         \    |       \      /
          '.   \       \   .'
            '-._|       \-'
                | |\     |
        __jgs___/ |_'.   /______''',

    r'''
   Waves are my home    .:~~--__                __--~~:.
 Wind is my life      ,:;'~'-,__~~--..,---..--~~__,-`~`::.
                    ,:;'        ''-,_ (. .)_,-``        `::.
                  ,;'                \ `\)/                `:.
                 '                    `--'                    `
   __._                       _.._                 _._
-~~    ~~--..__.._-~~~--..--~~    ~~--.__.---...-'~   ~~---...-.__seal__.''',

    r'''
                ________
            _.-'::'\____`.
          ,'::::'  |,------.
         /::::'    ||`-..___;
        ::::'      ||   / ___\
        |::       _||  [ [___]]
        |:   __,-'  `-._\__._/
        :_,-\  \| |,-'_,. . `.
        | \  \  | |.-'_,-\ \   ~
        | |`._`-| |,-|    \ \    ~
        |_|`----| ||_|     \ \     ~              _
        [_]     |_|[_]     [[_]      ~        __(  )
        | |    [[_]| |     `| |        ~    _(   )   )
        |_|    `| ||_|      |_|          ~ (    ) ) ))
        [_]     | |[_]      [_]          (_       _))
       /___\    [ ] __\    /___\           (( \   ) )
jrei          /___\                        (     ) )
                                             (  #  )''',

    r'''
              .==,_
             .===,_`\
           .====,_ ` \      .====,__
     ---     .==-,`~. \           `:`.__,
      ---      `~~=-.  \           /^^^   ...always on the go!
        ---       `~~=. \         /
                     `~. \       /
                       ~. \____./
              jgs        `.=====)
                      ___.--~~~--.__
            ___\.--~~~              ~~~---.._|/
            ~~~"                             /''',

    r'''
                   ~.
            Ya...___|__..ab.     .   .
             Y88b  \88b  \88b   (     )
              Y88b  :88b  :88b   `.oo'
              :888  |888  |888  ( (`-'
     .---.    d88P  ;88P  ;88P   `.`.
    / .-._)  d8P-"""|"""'-Y8P      `.`.
   ( (`._) .-.  .-. |.-.  .-.  .-.   ) )
    \ `---( O )( O )( O )( O )( O )-' /
     `.    `-'  `-'  `-'  `-'  `-'  .' CJ
       `---------------------------' ''',

    r'''
           .          .
 .          .                  .          .              .
       +.           _____  .        .        + .                    .
   .        .   ,-~"     "~-.                                +
              ,^ ___         ^. +                  .    .       .
             / .^   ^.         \         .      _ .
            Y  l  o  !          Y  .         __CL\H--.
    .       l_ `.___.'        _,[           L__/_\H' \\--_-          +
            |^~"-----------""~ ^|       +    __L_(=): ]-_ _-- -
  +       . !                   !     .     T__\ /H. //---- -       .
         .   \                 /               ~^-H--'
              ^.             .^            .      "       +.
                "-.._____.,-" .                    .
         +           .                .   +                       .
  +          .             +                                  .
         .             .      .       -Row
                                                        .''']

    assist_phrase = [
            ]

    find_headings = False # experimental
    find_tables = True
    max_math_length = 20

    if find_headings:
        html = u""
        for tag in soup.find("h4").next_siblings:
            if tag.name == "h4":
                continue
            else:
                html += str(tag)
        md = pypandoc.convert_text(html, 'plain', format='html')
        print(md)

    if find_tables:
        # find all tables
        for table in soup.findAll("table"):
            if table.findParent("table") is None:

                # remove long LaTeX strings
                for i in table.findAll('span', class_="math inline"):
                    if len(i.string)>max_math_length:
                        i.string = "(hidden math)" # we could also delete w. `i.decompose`

                # tag = text in upper left corner of table
                tag = table.thead.tr.findAll('th')[0].string
                if tag!=None:
                    md = pypandoc.convert_text(table, 'plain', format='html')
                    md = re.sub(r'\n+', '\n', md) # remove double newline
                    tables[tag] = md

        # store dictionary as JSON file 
        out = json.dumps(tables, indent=2)
        with open('tips.json', 'w') as f:
            f.write(out)

