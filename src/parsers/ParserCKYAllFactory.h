// -*- mode: c++ -*-
#ifndef _PARSERCKYALLFACTORY_H_
#define _PARSERCKYALLFACTORY_H_

#include "utils/ConfigTable.h"

#include "ParserCKYAllViterbi.h"
#include "ParserCKYAllMaxVar1B.h"
#include "ParserCKYAllMaxVarKB.h"
#include "ParserCKYAllMaxVarMultiple.h"

#include "utils/data_parsers/AnnotHistoriesParser.h"
#ifdef USE_THREADS
#include <tbb/task_scheduler_init.h>
#endif

namespace ParserCKYAllFactory {
    enum Parsing_Algorithm {MaxRule, Viterbi, MaxN, KMaxRule};
    ParserCKYAll * create_parser(ConfigTable& config);
    Parsing_Algorithm string_to_pa(const std::string& s);
};

#endif /* _PARSERCKYALLFACTORY_H_ */
