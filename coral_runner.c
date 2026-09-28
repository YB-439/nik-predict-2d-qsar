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
    if (argc < 2) {
        printf("Usage: coral_runner.exe <run_directory> [expected_smiles_count] [timeout_seconds]\n");
        return 1;
    }

    const char *run_dir = argv[1];
    int expected_count = (argc >= 3) ? atoi(argv[2]) : 1;
    int timeout_sec = (argc >= 4) ? atoi(argv[3]) : 30;

    char exe_path[MAX_PATH];
    char list_path[MAX_PATH];
    char out_path[MAX_PATH];

    snprintf(exe_path, sizeof(exe_path), "%s\\CORALSEA.exe", run_dir);
    snprintf(list_path, sizeof(list_path), "%s\\list.txt", run_dir);
    snprintf(out_path, sizeof(out_path), "%s\\ListModel.txt", run_dir);

    // Remove existing ListModel.txt before running
    DeleteFileA(out_path);

    // Launch CORALSEA.exe
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    if (!CreateProcessA(exe_path, NULL, NULL, NULL, FALSE, 0, NULL, run_dir, &si, &pi)) {
        fprintf(stderr, "Failed to launch %s (Error %lu)\n", exe_path, GetLastError());
        return 2;
    }

    // Wait for TForm1 window (up to 10 seconds)
    HWND target_hwnd = NULL;
    DWORD start_time = GetTickCount();
    while ((GetTickCount() - start_time) < 10000 && !target_hwnd) {
        FindWindowContext ctx = { pi.dwProcessId, NULL };
        EnumWindows(EnumWindowsCallback, (LPARAM)&ctx);
        if (ctx.hwnd) {
            target_hwnd = ctx.hwnd;
            break;
        }
        Sleep(200);
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
    while ((GetTickCount() - load_start) < 5000 && !btn_load) {
        btn_load = FindChildByPrefix(target_hwnd, "Load method");
        if (btn_load) break;
        Sleep(200);
    }

    if (btn_load) {
        PostMessageA(btn_load, BM_CLICK, 0, 0);
        Sleep(400);
    } else {
        fprintf(stderr, "Button 'Load method' not found\n");
    }

    // Click 'Import of current model'
    HWND btn_import = FindChildByPrefix(target_hwnd, "Import of current model");
    if (btn_import) {
        PostMessageA(btn_import, BM_CLICK, 0, 0);
        Sleep(500);
    } else {
        fprintf(stderr, "Button 'Import of current model' not found\n");
    }

    // Click 'Calculation model for a list of SMILES'
    HWND btn_calc = FindChildByPrefix(target_hwnd, "Calculation model for a list of SMILES");
    if (btn_calc) {
        PostMessageA(btn_calc, BM_CLICK, 0, 0);
    } else {
        fprintf(stderr, "Button 'Calculation model for a list of SMILES' not found\n");
    }

    // Wait for ListModel.txt to exist and contain results
    DWORD calc_start = GetTickCount();
    int success = 0;
    while ((GetTickCount() - calc_start) < (DWORD)(timeout_sec * 1000)) {
        Sleep(200);
        WIN32_FILE_ATTRIBUTE_DATA file_info;
        if (GetFileAttributesExA(out_path, GetFileExInfoStandard, &file_info)) {
            if (file_info.nFileSizeLow > 0) {
                success = 1;
                Sleep(200); // small delay to ensure flush
                break;
            }
        }
    }

    // Clean up process
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
