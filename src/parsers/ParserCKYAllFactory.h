// -*- mode: c++ -*-
#ifndef _PARSERCKYALLFACTORY_H_
#define _PARSERCKYALLFACTORY_H_

#include "utils/ConfigTable.h"
#include "ParserCKYAll.h"


namespace ParserCKYAllFactory {
    enum Parsing_Algorithm {MaxRule, Viterbi, MaxN, KMaxRule, MinDiv};
    ParserCKYAll * create_parser(ConfigTable& config);
    Parsing_Algorithm string_to_pa(const std::string& s);
};

#endif /* _PARSERCKYALLFACTORY_H_ */
