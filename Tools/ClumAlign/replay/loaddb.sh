#source ../analyzer/analyzer_install/env_set.sh
#export DB_DIR=/adaqfs/home/a-onl/siyu/DB/DB_PREXII_test
#export TREESEARCH=/adaqfs/home/a-onl/happexsp/TreeSearch
export DB_DIR=/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/ClumAlign/replay/DB
export MPDGEM=/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexAnalyzer


export LD_LIBRARY_PATH=$DB_DIR:$MPDGEM:$LD_LIBRARY_PATH
export PATH=$DB_DIR:$MPDGEM:$PATH
