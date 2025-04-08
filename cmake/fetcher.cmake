include(FetchContent)

function(git_fetch package version source required)
    if (version MATCHES "^[0-9]+\.([0-9]+\.)*[0-9]+$")
        if (${required})
            find_package(${package} ${version} REQUIRED)
        else()
            find_package(${package} ${version} QUIET)
        endif()
        # Add the version tag preceding v in case it wasn't found
        set(version "v${version}")
    endif()
    if (${${package}_FOUND})
        message(STATUS "Found system ${package}")
    else()
        message(STATUS "Fetch ${package} ${version} from ${source}")
        FetchContent_Declare(${package} GIT_REPOSITORY ${source} GIT_TAG ${version})
	FetchContent_MakeAvailable(${package})
        set(${package}_FOUND ON PARENT_SCOPE)
        set("${package}_SOURCE_DIR" "${${package}_SOURCE_DIR}" PARENT_SCOPE)
        set("${package}_BINARY_DIR" "${${package}_BINARY_DIR}" PARENT_SCOPE)
    endif()
endfunction()
