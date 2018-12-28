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
    tables['ascii'] = []

    tables['ascii'].append(r'''
    _                ___       _.--.
    \`.|\..----...-'`   `-._.-'_.-'`
    /  ' `         ,       __.--'
    )/' _/     \   `-_,   /
    `-'" `"\_  ,_.-;_.-\_ ',     fsc/as
        _.-'_./   {_.'   ; /
       {_.-``-'         {_/''')

    tables['ascii'].append(r'''
    __.-._
    '-._"7'
     /'.-c
     |  /T
snd _)_/LI''')

    tables['ascii'].append(r'''
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
                         \            \     `''')

    tables['ascii'].append(r'''
                  ,,__
        ..  ..   / o._)                   .---.
       /--'/--\  \-'||        .----.    .'     '.
      /        \_/ / |      .'      '..'         '-.
    .'\  \__\  __.'.'     .'          i-._
      )\ |  )\ |      _.'
     // \\ // \\
    ||_  \\|_  \\_
mrf '--' '--'' '--' ''')

    # http://www.oocities.org/spunk1111/wildlife.htm
    tables['ascii'].append(r'''
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
        __jgs___/ |_'.   /______''')

    tables['ascii'].append(r'''
   Waves are my home    .:~~--__                __--~~:.
 Wind is my life      ,:;'~'-,__~~--..,---..--~~__,-`~`::.
                    ,:;'        ''-,_ (. .)_,-``        `::.
                  ,;'                \ `\)/                `:.
                 '                    `--'                    `
   __._                       _.._                 _._
-~~    ~~--..__.._-~~~--..--~~    ~~--.__.---...-'~   ~~---...-.__seal__.''')

    tables['ascii'].append(r'''
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
                                             (  #  )''')

    tables['ascii'].append(r'''
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
            ~~~"                             /''')

    tables['ascii'].append(r'''
   ______________________________________________________________________
  | .     .               .   .                   .          .           |
  |              . .                     .      ___,_   _         .   .  |
  | .                       .      .          [:t_:::;t"t"+        .     |
  |      .                     .            . `=_ "`[ j.:\=\             |
  |             .      .              .        _,:-.| -"_:\=\  .         |
  |    .           .          .           _,-=":.:%.."+"+|:\=\        .  |
  |          .                   _ _____,:,,;,==.==+nnnpppppppt          |
  |                           _.;-^-._-:._::.'';nn;::m;:%%%%%%%\   .     |
  |  .       .              .;-'_::-:_"--;_:. ((888:(@) ,,;::^%%%,       |
  |                      __='::_:"`::::::::"-;_`YPP::; (d8B((@b."%\     .|
  |      ,------..    __,-:-:::::::::`::`::::::"--;_(@' 88P':^" ;nn:,    |
  |   ,-":%%%%::==.  ;-':::::`%%%\::---:::-:_::::::_"-;_.::((@,(88J::\   |
  |  /:::__ ::%::== """"""""""""""`------`.__.-:::::;___;;::`^__;;;:..7  |
  | /::.'  `.:%%=:=`-=,     . i                   .       """"           |
  |Y:::f    j :%%%%:::=::    ,^.    .        |-|  . .                    |
  |l   `.__+ :::%%%%:::_;[                        |o|                    |
  ||^~'-------------""~:^|                       _` ` _  .. __,,,,+++O#@@#
  |! ::::::::::%%%%==:{                       __j [,,j [#O|||O#@@#O++:|@@#
  | \ `::====: ==== :='            .__,,,++::::.j "  " [%+++::|@##O+::+O##
  |  \:== :: == :=='    __,,,+++|O|. +++..   :::j_[nnj_[_++:+%%_%%|+%|%+O#
  |   "-. =_:::: },+|O##|+::+|:++:::..    ::: .:+%%%%%%j [%O%%j [:+++|++|O
  |   _,,`-------' .+#O#+:||%+ ____   :: .. .:++|O###O%j `'  `' [:::::::++
  |.+:..:++|++||||+.O.::++:|::| _  |:...:++++|||+O##||%j [%..%j [+::LS:+%|
                              || \ |
                              ||  || eath Star II above Endor, surveyed by
                              ||_/ | the Super Star Destroyer "Executor".
                              |____|''')

    find_headings = False
    find_tables = True

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
        for table in soup.findAll("table"):
            if table.findParent("table") is None:
                tag = table.thead.tr.findAll('th')[0].string
                if tag!=None:
                    md = pypandoc.convert_text(table, 'plain', format='html')
                    md = re.sub(r'\n+', '\n', md)
                    tables[tag] = md

        out = json.dumps(tables, indent=2)
        with open('tips.json', 'w') as f:
            f.write(out)
