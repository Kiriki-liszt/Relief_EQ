cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

if(SMTG_MAC)
    set(vst3sdk_SOURCE_DIR "/Users/yoon/VST_SDK/vst3sdk")
elseif(SMTG_WIN)
    set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
endif(SMTG_MAC)

if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Relief_EQ
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.1.3.0 
    DESCRIPTION "Relief EQ VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Relief_EQ
    source/version.h
    source/RFEQ_cids.h
    source/RFEQ_svf.h           # added
    source/RFEQ_fft.h           # added
    source/RFEQ_fft.cpp         # added
    source/RFEQ_dataexchange.h  # added
    source/RFEQ_processor.h
    source/RFEQ_processor.cpp
    source/RFEQ_controller.h
    source/RFEQ_controller.cpp
    source/RFEQ_entry.cpp
)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(Relief_EQ
        PRIVATE
            resource/RFEQ_editor.uidesc
    )
    target_link_libraries(Relief_EQ
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Relief_EQ
        RESOURCES
            "resource/RFEQ_editor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (Relief_EQ
    RESOURCES
        resource/EB8EA8E4ADF058F5A94A0C179DFE6C8D_snapshot.png
        resource/EB8EA8E4ADF058F5A94A0C179DFE6C8D_snapshot_2.0x.png
)

target_link_libraries(Relief_EQ
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Relief_EQ)

if(SMTG_MAC)
    smtg_target_set_bundle(Relief_EQ
        BUNDLE_IDENTIFIER io.github.yg331.Relief.EQ
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(Relief_EQ
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Relief_EQ PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Relief_EQ)

        smtg_target_set_debug_executable(Relief_EQ
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)

# Add an AUv2 target
if (SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
    smtg_target_add_auv2(Relief_EQ-au
        BUNDLE_NAME Relief_EQ
        BUNDLE_IDENTIFIER io.github.yg331.Relief.EQ.audiounit
        INFO_PLIST_TEMPLATE resource/au-info.plist
        VST3_PLUGIN_TARGET Relief_EQ)
endif(SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
