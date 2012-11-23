#include "LorgApp.h"


LorgApp::LorgApp() : verbose(false), in(NULL), out(NULL)
#ifdef USE_THREADS
, tbb_task_scheduler(tbb::task_scheduler_init::deferred)
#endif
{}

LorgApp::~LorgApp()
{
  if (in != &std::cin) delete in;
  if (out != &std::cout) delete out;
}

bool LorgApp::init(int argc, char **argv)
{

  ConfigTable configuration(argc,argv,get_options());

  bool res = read_config(configuration);
  
#ifdef USE_THREADS
  unsigned nbthreads = configuration.get_value<unsigned>("nbthreads");
  if (nbthreads!=0)
    tbb_task_scheduler.initialize(nbthreads);
  else
    tbb_task_scheduler.initialize(tbb::task_scheduler_init::automatic);
#endif
  
  return res;
}


ConfigTable *  LorgApp::parse_config(int argc, char **argv)
{
  ConfigTable * configuration = new ConfigTable(argc,argv,get_options());

  if (!configuration) {
    std::cerr << "Unable to read options\n";
  }

  return configuration;
}


bool LorgApp::read_config(ConfigTable& configuration)
{
  verbose = configuration.exists("verbose");
  // parse config file if provided
  if(configuration.exists("config-file"))
    {
      if(verbose) std::clog << "Parsing configuration file." << std::endl;
      configuration.parse_config_file(configuration.get_value<std::string>("config-file"));
    }


  if (configuration.exists("help")) {
    configuration.print_help();
    return false;
  }

  return true;
}
