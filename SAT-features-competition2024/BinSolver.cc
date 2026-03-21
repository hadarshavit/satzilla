#include "BinSolver.h"

#include "global.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <signal.h>

BinSolver::BinSolver(const char *_name, int _argc, int _inputFileParam)
    : outFileCreated(false), argc(_argc),
      inputFileParam(_inputFileParam),
      name(_name)
{
  argv = new const char *[argc + 1];
}

BinSolver::~BinSolver()
{
  if (argv)
    delete[] argv;
}

/*
bool BinSolver::writeBinaryFile(char* nameBuf, const unsigned char* data, int length)
{
  int fd;

  sprintf(nameBuf, "%s", P_tmpdir);
  strcat(nameBuf, "/satzillaXXXXXX");
  nameBuf = mktemp(nameBuf);

  if ( !nameBuf ||
       (fd=open(nameBuf, O_CREAT | O_WRONLY | O_TRUNC, S_IRWXU))== -1 )
    {
      perror("Couldn't open solver for writing");
      return false;
    }



  if( write(fd, data, length)!= length )
    {
      perror("Writing solver file!");
      close(fd);
      unlink(nameBuf);
      return false;
    }

  close(fd);
  return true;
}
*/
int BinSolver::spawnBinary(const char *binFile, const char *const argv[], const char *outFileName, int timeout)
{
  // -- check timeout
  if (timeout == -1)
    timeout = (int)(gTimeOut - gSW.TotalLap());
  if (timeout < 1) // -- no point in running
    return 1;
  if (DEB)
    printf("c spawning binary file with %d seconds ...\n", timeout);

  int pid = fork();

  if (pid == -1)
  {
    perror("Forking Child!");
    return 123;
  }

  if (pid)
  {
    // -- parent
    int status;
    wait(&status);

    if (WIFSIGNALED(status))
    {
      printf("c child exit by a signal\n");
      return 128 + WTERMSIG(status);
    }
    else if (WIFEXITED(status))
    {
      if (DEB)
        printf("c child exited successfully\n");
      return WEXITSTATUS(status);
    }
    // -- else
    printf("c child quit for some reason\n");
    return 123;
  }
  else
  {
    // -- child
    //-- try to do redirection
    int fd = open(outFileName, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
    if (fd == -1)
    {
      perror("Redirect failed, results undefined!");
    }
    else
    {
      dup2(fd, 1);
      close(fd);
    }

    // -- go till hard kill, just in case
    signal(SIGXCPU, SIG_IGN);

    // -- setup timeout
    struct rlimit rlim;
    rlim.rlim_cur = rlim.rlim_max = timeout;

    int res = setrlimit(RLIMIT_CPU, &rlim);

    printf("c child timeout set to %d\n", timeout);

    // -- finally, execute
    char **execArgv = new char *[argc + 1];
    for (int i = 0; i < argc; ++i)
      execArgv[i] = const_cast<char *>(argv[i]);
    execArgv[argc] = nullptr;
    execv(const_cast<char *>(binFile), execArgv);
    perror("Exec");
    exit(123);
  }

  return 123;
}

int BinSolver::execute(const char *inputFile, int timeout)
{
  if (DEB)
    printf("c Bin: Executing %s\n", name);
  /* -- First, spit out executable
  solverFileCreated = writeBinaryFile(solverName, data, dataSize);
  if(!solverFileCreated)
    return -1;
  */
  // -- Create outfile
  snprintf(outFileName, sizeof(outFileName), "%s/outputXXXXXX", P_tmpdir);
  int fd = mkstemp(outFileName);
  if (fd == -1)
  {
    perror("mkstemp failed for solver output file");
    return 123;
  }
  close(fd);
  outFileCreated = true;

  // -- execute!
  snprintf(solverName, sizeof(solverName), "%s/satzilla_Solvers/%s", mypath, name);
  // printf("c Bin: Executing %s\n", solverName);
  argv[inputFileParam] = inputFile;
  argv[0] = solverName;
  int res = spawnBinary(solverName, argv, outFileName, timeout);

  // -- try to cleanup

  return res;
}

void BinSolver::cleanup()
{
  if (outFileCreated)
  {
    if (unlink(outFileName) == 0)
      outFileCreated = false;
  }
}
