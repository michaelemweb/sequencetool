SET(SOURCES
  Main.C
  Utils.C
  NtAlign.C
  FastaConvert.C
  ConsensusAlign.C
  MakeConsensus.C
  ssw_cpp.cpp
  ssw.c
)

ADD_EXECUTABLE(sequencetool ${SOURCES})

TARGET_LINK_LIBRARIES(sequencetool seq)

INSTALL_TARGETS(/bin sequencetool)
