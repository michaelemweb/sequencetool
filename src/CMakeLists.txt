SET(SOURCES
  Mpf.C
  Utils.C
  NtAlign.C
  FastaConvert.C
)

INCLUDE_DIRECTORIES(/opt/software/boost1_43/include)

ADD_EXECUTABLE(mpf-sequencetool ${SOURCES})

TARGET_LINK_LIBRARIES(mpf-sequencetool seq)

INSTALL_TARGETS(/bin mpf-sequencetool)
