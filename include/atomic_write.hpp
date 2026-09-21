#ifndef ATOMIC_WRITE_HPP
#define ATOMIC_WRITE_HPP

#include <cstdio>
#include <iostream>
#include <string>

// Move a finished temporary file onto its final name.
//
// The outputs here are read by a separate process while the simulation runs,
// and it polls for them by name. Writing straight to the final path publishes
// the name as soon as the file is created, long before the contents are
// complete, so the reader can pick up a partial file that still looks
// structurally valid. Writing to a temporary name and renaming afterwards
// closes that window: rename(2) is atomic within a directory, so the final name
// only ever resolves to a finished file.
//
// The rename is for visibility to a concurrent reader, which sees the same page
// cache, so there is no fsync here; this does not make the write durable across
// a machine crash.
inline int publish_temp_file(const std::string &temp_filename,
			     const std::string &final_filename,
			     const std::string &context)
{
  if (std::rename(temp_filename.c_str(),final_filename.c_str()) != 0)
    {
      std::cout << "ERROR: could not rename " << temp_filename
		<< " to " << final_filename
		<< " in " << context << std::endl;
      std::remove(temp_filename.c_str());
      return 1;
    }

  return 0;
}

#endif
