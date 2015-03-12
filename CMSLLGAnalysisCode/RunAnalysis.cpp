#include "LLGAnalysis.h"
#include <vector>
#include <string>

int main( int argc, char **argv ) {

    vector<string> fileNames;
    for( int iFile = 1; iFile < argc; ++iFile ) {
        fileNames.push_back( string(argv[iFile]) );
    }
    LLGAnalysis *analysis = LLGAnalysis::GetInstance(fileNames, "RecoData" );
    analysis->Init();
    analysis->RunEventLoop();
    analysis->FinishRun();
    return 0;
}
