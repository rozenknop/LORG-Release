add_subdirectory(training)

set(common_sources 
         rules/Production.cpp
         rules/Rule.cpp
         grammars/compact_binary_rules.cpp
         grammars/Grammar.cpp
         utils/SymbolTable.cpp
         utils/Tagger.cpp
         utils/TokeniserFactory.cpp
         utils/Tokeniser.cpp
         utils/PtbPsTree.cpp
         utils/LorgConstants.h
         utils/ConfigTable.cpp
         utils/RandomGenerator.cpp
         lexicon/WordSignatureFactory.cpp
         lexicon/BerkeleyEnglishUMapping.cpp
         lexicon/BaselineFrenchUMapping.cpp
         lexicon/EnglishIGMapping.cpp
         lexicon/FrenchIGMapping.cpp
         lexicon/ItalianIGMapping.cpp
         lexicon/ArabicUMapping.cpp
         lexicon/ArabicIGMapping.cpp
         Word.cpp
         LorgApp.cpp
         LorgParseApp.cpp
         feature_extract/Extract.cpp
         feature_extract/Feature.cpp
         ParseSolution.cpp
   )

set(simple_lorgparser_SOURCES
         ${common_sources}
         simple_lorgparser.cpp
         grammars/compact_binary_rules_simple.cpp
         utils/data_parsers/RuleInputParser.cpp
         edges/Edge.cpp
         parsers/ParserCKYBest.cpp
         PCKYBestCell.cpp
         SimpleLorgParseApp.cpp)

set(twostage_lorgparser_SOURCES
         ${common_sources}
         twostage_lorgparser.cpp
         TwoStageLorgParseApp.cpp
         edges/ViterbiProbability.cpp
         edges/MaxRuleProbability1B.cpp
         edges/MaxRuleProbabilityKB.cpp
         edges/MaxRuleProbabilityMultiple.cpp
         parsers/ParserCKYAll.cpp
         parsers/ParserCKYAllFactory.cpp
         parsers/ParserCKYAllViterbi.cpp
         parsers/ParserCKYAllMaxVar1B.cpp
         parsers/ParserCKYAllMaxVarKB.cpp
         parsers/ParserCKYAllMaxVarMultiple.cpp
         parsers/ParserCKYAllMinDivKB.cpp
         grammars/compact_binary_rules_annotated.cpp
         grammars/AnnotatedLabelsInfo.cpp
         grammars/GrammarAnnotated.cpp
         utils/data_parsers/BURuleInputParser.cpp
         utils/data_parsers/AnnotHistoriesParser.cpp
         edges/AnnotationInfo.cpp
         rules/BRule.cpp rules/URule.cpp rules/LexicalRule.cpp)

set(common_include_dir ${CMAKE_CURRENT_SOURCE_DIR} ${TBB_INCLUDE_DIRS})
set(common_libs  ${TBB_LIBS} ${RTLIB} ${TCMALLOC_LIB})

#message(warning "In src : Include directories are : ${common_include_dir}")
# message(warning "In src : Common libs are : ${common_libs}")

include_directories(${common_include_dir})

add_executable(simple_lorgparser ${simple_lorgparser_SOURCES})
add_executable(twostage_lorgparser ${twostage_lorgparser_SOURCES})

target_link_libraries(simple_lorgparser ${common_libs} ${Boost_REGEX_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(twostage_lorgparser ${common_libs} ${Boost_LIBRARIES})



# install
install (TARGETS simple_lorgparser twostage_lorgparser DESTINATION bin)
