// -*- mode: c++ -*-
#ifndef PARSESOLUTION_H
#define PARSESOLUTION_H

#include "utils/PtbPsTree.h"
#include "feature_extract/Extract.h"


/////////////////////////////////////
template <typename IdentifierType,
          class IgnoredAbstractProduct>
class DefaultFactoryError
{
public:
  class Exception : std::exception
  {
  public:
    Exception(const IdentifierType& id) : unknown_id(id) {};

    const IdentifierType& get_id() {return unknown_id;};

    virtual const char * what() { return "Unknown Object Type passed to Factory";};

    ~Exception() throw() {};


  private:
    IdentifierType unknown_id;
  };

protected:
  IgnoredAbstractProduct * on_unknown_type_id(const IdentifierType&id)
  {
    throw Exception(id);
  };
};


template <class AbstractProduct,
          typename IdentifierType,
          typename DataIn,
          typename ProductCreator = AbstractProduct* (*)(const DataIn&),
          template <typename, class> class FactoryErrorPolicy = DefaultFactoryError
          >
class Factory : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
{
public:
  bool register_type(const IdentifierType& id, ProductCreator creator)
   {
    return assoc.insert(typename AssocMap::value_type(id,creator)).second;
   };

  bool unregister_type(const IdentifierType& id)
  {
    return assoc.erase(id) == 1;
  };

  AbstractProduct* create_object(const IdentifierType& id, const DataIn& data)
  {
    typename AssocMap::const_iterator i = assoc.find(id);
    if(i != assoc.end())
      {
        return (i->second)(data);
      }
    else
      {
        return this->on_unknown_type_id(id);
      }
  };

 private:
   typedef std::map<IdentifierType, ProductCreator> AssocMap;
   AssocMap assoc;
 };

////////////////////////////////////

class ParseSolutionFactory;

struct parse_solution
{
  enum parse_solution_format { UNIX, JSON};


  const std::string& sentence;
  int id_sentence;
  unsigned length;
  std::vector<std::pair<PtbPsTree *, double> >& trees;
  double time;
  bool verbose;

  const std::vector<std::string>& comments;
  bool extract_features;


  static parse_solution_format format_from_string(const std::string&);

  static Extract extractor;
  static void init_feature_extractor();
  static ParseSolutionFactory factory;

  virtual ~parse_solution() {};


  virtual std::ostream& print(std::ostream& out) const {return out;};


  parse_solution(const std::string& s, int id, unsigned l, std::vector<std::pair<PtbPsTree*, double> >& trs, const double& ti, bool ve, const std::vector<std::string>& co, bool ef)
    : sentence(s), id_sentence(id), length(l), trees(trs), time(ti), verbose(ve), comments(co), extract_features(ef)
  {};

  parse_solution(const parse_solution& other)
    : sentence(other.sentence),
      id_sentence(other.id_sentence),
      length(other.length),
      trees(other.trees),
      time(other.time),
      verbose(other.verbose),
      comments(other.comments),
      extract_features(other.extract_features)
  {};
};


class ParseSolutionFactory : public Factory<parse_solution,int,parse_solution>
{
};

struct unix_parse_solution : parse_solution
{
  unix_parse_solution(const parse_solution& p) : parse_solution(p) {};

  static parse_solution * create(const parse_solution&p) {return new unix_parse_solution(p);}
  static bool init() {return parse_solution::factory.register_type(UNIX, create);}
  std::ostream& print(std::ostream& out) const;
};
std::ostream& operator<<(std::ostream& out, const unix_parse_solution& ps);


struct json_parse_solution : parse_solution
{
  json_parse_solution(const parse_solution& p) : parse_solution(p) {};

  static parse_solution * create(const parse_solution&p) {return new json_parse_solution(p);}
  static bool init() {return parse_solution::factory.register_type(JSON, create);}
  std::ostream& print(std::ostream& out) const;
};

#endif //PARSESOLUTION_H
