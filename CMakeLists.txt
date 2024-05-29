cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(TestEQ_01
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "TestEQ_01 VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(TestEQ_01
    source/version.h
    source/RFEQ_cids.h
    source/RFEQ_svf.h           # added
    source/RFEQ_dataexchange.h  # added
    source/RFEQ_processor.h
    source/RFEQ_processor.cpp
    source/RFEQ_controller.h
    source/RFEQ_controller.cpp
    source/RFEQ_entry.cpp
)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(TestEQ_01
        PRIVATE
            resource/RFEQ_editor.uidesc
    )
    target_link_libraries(TestEQ_01
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(TestEQ_01
        RESOURCES
            "resource/RFEQ_editor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (TestEQ_01
    RESOURCES
        resource/EB8EA8E4ADF058F5A94A0C179DFE6C8D_snapshot.png
        resource/EB8EA8E4ADF058F5A94A0C179DFE6C8D_snapshot_2.0x.png
)

target_link_libraries(TestEQ_01
    PRIVATE
        sdk
)

smtg_target_configure_version_file(TestEQ_01)

if(SMTG_MAC)
    smtg_target_set_bundle(TestEQ_01
        BUNDLE_IDENTIFIER io.github.yg331.RFEQ
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(TestEQ_01
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(TestEQ_01 PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT TestEQ_01)

        smtg_target_set_debug_executable(TestEQ_01
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)
