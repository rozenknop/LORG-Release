// -*- mode: c++ -*-
#ifndef TWOSTAGELORGPARSEAPP_H
#define TWOSTAGELORGPARSEAPP_H

#include "LorgParseApp.h"
#include "ParseSolution.h"

class ParserCKYAll;

class TwoStageLorgParseApp : public LorgParseApp
{
public:
  TwoStageLorgParseApp();
  ~TwoStageLorgParseApp();
  int run();

private:
  bool read_config(ConfigTable& configuration);
  LorgOptions get_options() const;

  ParserCKYAll * parser;
  Tagger tagger;
  bool output_annotations;
  unsigned kbest;

  bool extract_features;
  parse_solution::parse_solution_format output_format;
};

#endif // TWOSTAGELORGPARSEAPP_H
