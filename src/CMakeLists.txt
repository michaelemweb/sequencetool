SET(SOURCES
  Mpf.C
  Utils.C
  NtAlign.C
  FastaConvert.C
)

#INCLUDE_DIRECTORIES(/opt/software/include/boost-1_37)

ADD_EXECUTABLE(mpf-sequencetool ${SOURCES})

TARGET_LINK_LIBRARIES(mpf-sequencetool seq)

INSTALL_TARGETS(/bin mpf-sequencetool)
