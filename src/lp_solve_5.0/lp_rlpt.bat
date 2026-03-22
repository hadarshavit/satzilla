flex -L -l lp_rlpt.l
sed -e "s/yy/lpt_yy/g" lex.yy.c >lp_rlpt.h
del lex.yy.c

bison --no-lines -y lp_rlpt.y
sed -e "s/yy/lpt_yy/g" y.tab.c >lp_rlpt.c
del y.tab.c
