set nocompatible
autocmd BufReadPost * if line("'\"") && line("'\"") <= line("$") | exe "normal `\"" | endif
colorscheme torte
"colorscheme elflord
"colorscheme evening
"colorscheme yeller
"set termcap

set expandtab
set shiftwidth=2
set tabstop=2
set bs=2                " allow backspacing over everything in insert mode set bs=2
"set cursorline          " highlight cursor line
set nonu                  " set line number
set hlsearch            " highlight search

"ses guitablabel=%{tabpagenr()}.%t\ %m
"set spell
set mouse=""

set viminfo='20,\"50    " read/write a .viminfo file, don't store more
                        " than 50 lines of registers
set history=100          " keep 50 lines of command line history
set ruler               " show the cursor position all the time
set backupdir=/afs/ihep.ac.cn/users/w/wenlj/backup
syntax on

set sm
set ai                  " always set autoindenting on set ai
set cindent             " auto indent
set cinoptions={0,1s,t0,n-2,p2s,(03s,=.5s,>1s,=1s,:1s
set t_ti=[?47h        " restore screen after quit
set t_te=[?47l
"au BufEnter,BufRead *.py setlocal smartindent cinwords=if,elif,else,for,while,try,except,finally,def,class

"if &term == "xterm"
  set t_kb=^?
  set t_Co=8
  set t_Sb=^[[4%dm
  set t_Sf=^[[3%dm
"endif
