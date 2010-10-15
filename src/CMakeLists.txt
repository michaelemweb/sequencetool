SET(SOURCES
  Mpf.C
  NtAlign.C
)

#INCLUDE_DIRECTORIES(/opt/software/include/boost-1_37)

ADD_EXECUTABLE(mpf-sequencetool ${SOURCES})

TARGET_LINK_LIBRARIES(mpf-sequencetool seq boost_regex-mt)

INSTALL_TARGETS(/bin mpf-sequencetool)
