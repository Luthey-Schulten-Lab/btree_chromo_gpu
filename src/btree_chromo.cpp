#include <btree_chromo.hpp>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include <fcntl.h>
#include <unistd.h>
#include "fused_cg_minimize.h"

// joins the CUDA context pre-warm thread on every return from main
struct wcm_prewarm_guard { ~wcm_prewarm_guard() { fused_prewarm_join(); } };

// btree_chromo --serve runs one directive file per stdin line, each through a fresh btree_driver exactly as a separate
// btree_chromo process would, so the per-hook process start, CUDA context, MPI and Kokkos initialisation are paid once per run.
// The regular output goes to /dev/null (the 4DWCM launches btree with stdout to DEVNULL; LAMMPS keeps its per-hook -log file);
// one protocol line per request goes to the original stdout: "BTREE_SERVE_DONE <rc> <seconds> <directive file>". EXIT or EOF ends.
extern bool wcm_serve_mode;
void wcm_serve_finalize();
static int wcm_run_directives(const std::string &drctv_filename);

static int wcm_serve()
{
  wcm_serve_mode = true;
  const int proto = dup(1);
  const int devnull = open("/dev/null", O_WRONLY);
  if (proto < 0 || devnull < 0) return 1;
  dup2(devnull, 1);
  close(devnull);
  FILE *pf = fdopen(proto, "w");
  fprintf(pf, "BTREE_SERVE_READY\n");
  fflush(pf);
  std::string line;
  while (std::getline(std::cin, line))
    {
      if (line.empty()) continue;
      if (line == "EXIT") break;
      // optional "<directive file>\t<log file>" (tests): that request's stdout and stderr go to the file, like a replay's 2>&1
      std::string path = line, log;
      const size_t tab = line.find('\t');
      if (tab != std::string::npos) { path = line.substr(0, tab); log = line.substr(tab + 1); }
      int saved = -1, saved2 = -1;
      if (!log.empty())
        {
          const int lf = open(log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
          if (lf >= 0) { std::cout.flush(); std::cerr.flush(); fflush(stdout); fflush(stderr); saved = dup(1); saved2 = dup(2); dup2(lf, 1); dup2(lf, 2); close(lf); }
        }
      const auto t0 = std::chrono::steady_clock::now();
      int rc;
      try { rc = wcm_run_directives(path); }
      catch (const std::exception &e) { std::cerr << "btree_chromo --serve exception: " << e.what() << std::endl; rc = 2; }
      catch (...) { rc = 3; }
      std::cout.flush();
      fflush(stdout);
      if (saved >= 0) { std::cerr.flush(); fflush(stderr); dup2(saved, 1); close(saved); dup2(saved2, 2); close(saved2); }
      const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      fprintf(pf, "BTREE_SERVE_DONE %d %.3f %s\n", rc, sec, path.c_str());
      fflush(pf);
    }
  wcm_serve_finalize();
  return 0;
}

int main(int argc, char *argv[])
{
  fused_prewarm_cuda_context_async();
  wcm_prewarm_guard wcm_prewarm;

  if (argc != 2)
    {
      std::cout << "ERROR: missing directive file as command-line argument" << std::endl;
      return 1;
    }

  if (std::string(argv[1]) == "--serve") return wcm_serve();
  return wcm_run_directives(std::string(argv[argc-1]));
}

static int wcm_run_directives(const std::string &drctv_filename)
{
  int error_code = 0;
  btree_driver driver;

  // read the directives from the directive file
  driver.read_directives(drctv_filename);

  // print the directives to be executed
  driver.print_directives();

  // parse the directives into commands and parameters
  driver.parse_directives();

  // validate the numbers of parameters
  error_code = driver.validate_command_sequence_parameters();
  if (error_code != 0) return 1;
  
  // expand the metacommands
  error_code = driver.expand_metacommands();
  if (error_code != 0) return 1;

  // validate the command sequence
  error_code = driver.validate_command_sequence();
  if (error_code != 0) return 1;

  // print the set of commands
  driver.print_commands();
    
  // execute the commands
  error_code = driver.execute_commands();
  if (error_code != 0)
    {
      std::cout << "error during command execution" << std::endl;
      return 1;
    } 

  return 0;
  
}
