SET(SOURCES
  Mpf.C
  Utils.C
  NtAlign.C
  FastaConvert.C
  ConsensusAlign.C
  MakeConsensus.C
  ssw_cpp.cpp
  ssw.c
)

ADD_EXECUTABLE(mpf-sequencetool ${SOURCES})

TARGET_LINK_LIBRARIES(mpf-sequencetool seq)

INSTALL_TARGETS(/bin mpf-sequencetool)
