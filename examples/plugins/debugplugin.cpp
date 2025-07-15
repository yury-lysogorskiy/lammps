
#include "lammpsplugin.h"

#include "comm.h"
#include "command.h"
#include "error.h"
#include "version.h"

#include <cfenv>
#include <cstring>

namespace LAMMPS_NS {
class DebugFP : public Command {
 public:
  DebugFP(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **) override;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;

void DebugFP::command(int argc, char **argv)
{
  if (argc < 1) utils::missing_cmd_args(FLERR, "debugfp", error);

  if (strcmp(argv[0],"default") == 0) {
    fesetenv(FE_DFL_ENV);
    return;
  } else if (strcmp(argv[0],"trap") == 0) {
    fedisableexcept(FE_ALL_EXCEPT);
    for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "divbyzero") == 0) {
        feenableexcept(FE_DIVBYZERO);
      } else if (strcmp(argv[i], "inexact") == 0) {
        feenableexcept(FE_INEXACT);
      } else if (strcmp(argv[i], "invalid") == 0) {
        feenableexcept(FE_INVALID);
      } else if (strcmp(argv[i], "overflow") == 0) {
        feenableexcept(FE_OVERFLOW);
      } else if (strcmp(argv[i], "underflow") == 0) {
        feenableexcept(FE_UNDERFLOW);
      } else {
        error->all(FLERR, i, "Unknown floating-point trap {}", argv[i]);
        return;
      }
    }
    return;
  } else {
   error->all(FLERR, Error::ARGZERO, "Unknown debugfp keyword {}", argv[0]);
  }
}

static Command *debugcreator(LAMMPS *lmp)
{
  return new DebugFP(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "command";
  plugin.name = "debugfp";
  plugin.info = "Debug Floating-Point command v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &debugcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
