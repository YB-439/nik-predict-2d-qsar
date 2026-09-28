/*
 * coral_runner.c - Headless Win32 driver for CORALSEA.exe
 * Compiles with: i686-w64-mingw32-gcc -O2 coral_runner.c -o coral_runner.exe
 * Runs natively on Windows or under Wine on Linux (Xvfb)
 */

#include <windows.h>
#include <stdio.h>
#include <string.h>

#define BM_CLICK 0x00F5

typedef struct {
    DWORD pid;
    HWND hwnd;
} FindWindowContext;

static BOOL CALLBACK EnumWindowsCallback(HWND hwnd, LPARAM lParam) {
    FindWindowContext *ctx = (FindWindowContext *)lParam;
    DWORD pid = 0;
    GetWindowThreadProcessId(hwnd, &pid);
    if (pid == ctx->pid) {
        char cls[256];
        if (GetClassNameA(hwnd, cls, sizeof(cls))) {
            if (strcmp(cls, "TForm1") == 0) {
                ctx->hwnd = hwnd;
                return FALSE; // Stop enumeration
            }
        }
    }
    return TRUE;
}

typedef struct {
    const char *prefix;
    HWND hwnd;
} FindChildContext;

static BOOL CALLBACK EnumChildCallback(HWND hwnd, LPARAM lParam) {
    FindChildContext *ctx = (FindChildContext *)lParam;
    char text[256];
    if (GetWindowTextA(hwnd, text, sizeof(text))) {
        if (strncmp(text, ctx->prefix, strlen(ctx->prefix)) == 0) {
            ctx->hwnd = hwnd;
            return FALSE; // Found
        }
    }
    return TRUE;
}

static HWND FindChildByPrefix(HWND parent, const char *prefix) {
    FindChildContext ctx = { prefix, NULL };
    EnumChildWindows(parent, EnumChildCallback, (LPARAM)&ctx);
    return ctx.hwnd;
}

int main(int argc, char *argv[]) {
    const char *run_dir = (argc >= 2) ? argv[1] : ".";
    int timeout_sec = (argc >= 4) ? atoi(argv[3]) : 25;

    // Change to target directory if specified
    if (run_dir && strcmp(run_dir, ".") != 0) {
        SetCurrentDirectoryA(run_dir);
    }

    const char *exe_path = "CORALSEA.exe";
    const char *out_path = "ListModel.txt";

    // Remove existing ListModel.txt before running
    DeleteFileA(out_path);

    // Launch CORALSEA.exe in current directory
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    if (!CreateProcessA(NULL, (LPSTR)"CORALSEA.exe", NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
        fprintf(stderr, "Failed to launch CORALSEA.exe (Error %lu)\n", GetLastError());
        return 2;
    }

    // Wait for TForm1 window (up to 8 seconds)
    HWND target_hwnd = NULL;
    DWORD start_time = GetTickCount();
    while ((GetTickCount() - start_time) < 8000 && !target_hwnd) {
        FindWindowContext ctx = { pi.dwProcessId, NULL };
        EnumWindows(EnumWindowsCallback, (LPARAM)&ctx);
        if (ctx.hwnd) {
            target_hwnd = ctx.hwnd;
            break;
        }
        Sleep(150);
    }

    if (!target_hwnd) {
        fprintf(stderr, "Timed out waiting for CORALSEA TForm1 window\n");
        TerminateProcess(pi.hProcess, 1);
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
        return 3;
    }

    // Wait for and click 'Load method'
    HWND btn_load = NULL;
    DWORD load_start = GetTickCount();
    while ((GetTickCount() - load_start) < 4000 && !btn_load) {
        btn_load = FindChildByPrefix(target_hwnd, "Load method");
        if (btn_load) break;
        Sleep(150);
    }

    if (btn_load) {
        PostMessageA(btn_load, BM_CLICK, 0, 0);
        Sleep(350);
    }

    // Click 'Import of current model'
    HWND btn_import = FindChildByPrefix(target_hwnd, "Import of current model");
    if (btn_import) {
        PostMessageA(btn_import, BM_CLICK, 0, 0);
        Sleep(400);
    }

    // Click 'Calculation model for a list of SMILES'
    HWND btn_calc = FindChildByPrefix(target_hwnd, "Calculation model for a list of SMILES");
    if (btn_calc) {
        PostMessageA(btn_calc, BM_CLICK, 0, 0);
    }

    // Wait for ListModel.txt to exist and contain results
    DWORD calc_start = GetTickCount();
    int success = 0;
    while ((GetTickCount() - calc_start) < (DWORD)(timeout_sec * 1000)) {
        Sleep(150);
        WIN32_FILE_ATTRIBUTE_DATA file_info;
        if (GetFileAttributesExA(out_path, GetFileExInfoStandard, &file_info)) {
            if (file_info.nFileSizeLow > 0) {
                success = 1;
                Sleep(150);
                break;
            }
        }
    }

    // Clean up process immediately
    TerminateProcess(pi.hProcess, 0);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    if (success) {
        printf("SUCCESS: ListModel.txt generated successfully.\n");
        return 0;
    } else {
        fprintf(stderr, "Timed out waiting for ListModel.txt\n");
        return 4;
    }
}
