import os
import sys
import socket
import webbrowser
import threading
import time
import subprocess
import ctypes
from ctypes import wintypes
import uvicorn

# Fix stdout and stderr encoding on Windows
if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if sys.stderr and hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# Set console title & disable QuickEdit mode (to prevent freezing when user clicks in console)
if os.name == "nt":
    try:
        ctypes.windll.kernel32.SetConsoleTitleW("NIK-Predict: 2D-QSAR Local Host Server (Port 8005)")
    except Exception:
        pass
    try:
        kernel32 = ctypes.windll.kernel32
        hStdin = kernel32.GetStdHandle(-10)  # STD_INPUT_HANDLE
        mode = wintypes.DWORD()
        if kernel32.GetConsoleMode(hStdin, ctypes.byref(mode)):
            # 0x0040 is ENABLE_QUICK_EDIT_MODE, 0x0080 is ENABLE_EXTENDED_FLAGS
            kernel32.SetConsoleMode(hStdin, (mode.value & ~0x0040) | 0x0080)
    except Exception:
        pass

# Windows Anti-Sleep Constants
ES_CONTINUOUS = 0x80000000
ES_SYSTEM_REQUIRED = 0x00000001
ES_AWAYMODE_REQUIRED = 0x00000040

def activate_anti_sleep_mode() -> bool:
    """
    Prevents the host computer from entering sleep or standby mode while the server is active.
    This ensures that even after 1 hour or overnight, the local host continues running without suspension.
    """
    if os.name != "nt":
        return False
    try:
        kernel32 = ctypes.windll.kernel32
        flags = ES_CONTINUOUS | ES_SYSTEM_REQUIRED | ES_AWAYMODE_REQUIRED
        kernel32.SetThreadExecutionState(flags)

        def _heartbeat():
            while True:
                time.sleep(30)
                try:
                    kernel32.SetThreadExecutionState(flags)
                except Exception:
                    pass

        t = threading.Thread(target=_heartbeat, daemon=True)
        t.start()
        return True
    except Exception:
        return False

def check_vcredist_installed() -> bool:
    """
    Checks if Microsoft Visual C++ 2015-2022 Redistributable (x64) is installed.
    RDKit and C++ scientific binaries depend on this.
    """
    # If rdkit already imports cleanly, runtime is guaranteed functional
    try:
        import rdkit
        return True
    except Exception:
        pass

    if os.name != "nt":
        return True

    # 1. Check registry
    try:
        import winreg
        keys = [
            r"SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\X64",
            r"SOFTWARE\WOW6432Node\Microsoft\VisualStudio\14.0\VC\Runtimes\X64",
        ]
        for k in keys:
            try:
                with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, k) as reg_key:
                    val, _ = winreg.QueryValueEx(reg_key, "Installed")
                    if val == 1:
                        return True
            except OSError:
                pass
    except Exception:
        pass

    # 2. Check DLLs in System32
    sys32 = os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32")
    if os.path.isfile(os.path.join(sys32, "vcruntime140.dll")) and os.path.isfile(os.path.join(sys32, "msvcp140.dll")):
        return True

    return False

def find_available_port(preferred_port=8005):
    candidates = [preferred_port, 8006, 8007, 8008, 8080, 8888]
    for p in candidates:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("0.0.0.0", p))
                return p
            except OSError:
                continue
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("0.0.0.0", 0))
        return s.getsockname()[1]

def run_system_diagnostics():
    """
    Performs startup verification checks in the console window.
    """
    print("\n  [SYSTEM DIAGNOSTICS & ENVIRONMENT SELF-CHECK]")
    
    # 1. VC++ Redistributable check
    if check_vcredist_installed():
        print("  [OK] Microsoft Visual C++ Runtime: Verified (x64 detected)")
    else:
        print("  [!] Notice: Visual C++ Runtime not detected in standard paths.")

    # 2. RDKit cheminformatics engine check
    try:
        import rdkit
        from rdkit import Chem
        print("  [OK] RDKit Cheminformatics Engine: Verified")
    except Exception as e:
        print(f"  [!] Notice: Cheminformatics engine import: {e}")

    # 3. Bundled models check
    base_dir = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    models_dir = os.path.join(base_dir, "models")
    if os.path.isdir(models_dir):
        r1 = os.path.isdir(os.path.join(models_dir, "Run-1"))
        r2 = os.path.isdir(os.path.join(models_dir, "Run-2"))
        r3 = os.path.isdir(os.path.join(models_dir, "Run-3"))
        if r1 and r2 and r3:
            print("  [OK] CORALSEA Models: Verified (Run-1, Run-2, Run-3 bundled)")
        else:
            print("  [!] Some CORALSEA model subfolders are missing.")
    else:
        print("  [!] Models directory not detected in bundle.")

    # 4. Anti-sleep activation
    if activate_anti_sleep_mode():
        print("  [OK] Host Anti-Sleep Guard: ACTIVE (Prevents standby; continuous host mode)")
    else:
        print("  [!] Anti-sleep mode not available on this platform.")

def main():
    try:
        print("=" * 76)
        print("  Multi-Kinase 2D-QSAR Prediction Platform (CORALSEA Monte Carlo)")
        print("  Drug Design Synthesis Lab")
        print("  Department of Pharmaceutical Sciences and Drug Research")
        print("  Punjabi University, Patiala, Punjab, India")
        print("  Contact: drugdesignsynthesislab@gmail.com")
        print("=" * 76)

        # Run startup diagnostics & dependency check
        run_system_diagnostics()

        env_port = os.environ.get("PORT")
        if env_port:
            port = int(env_port)
        else:
            port = find_available_port(8005)

        host = os.environ.get("HOST", "0.0.0.0")
        display_url = f"http://localhost:{port}/"

        print(f"\n  [OK] Network Port Assigned: {port}")
        print(f"  Starting Web Portal at: {display_url}")
        print("  Active Kinase Engine: NF-kB Inducing Kinase (NIK / MAP3K14)")
        print("  Browser will open automatically.")
        print("  Keep this CMD window open while using the tool.")
        print("  Press Ctrl+C to terminate the server.\n" + "-" * 76)

        # Launch browser automatically in separate daemon thread
        def _open_browser():
            time.sleep(1.0)
            try:
                webbrowser.open(display_url)
            except Exception:
                pass

        threading.Thread(target=_open_browser, daemon=True).start()

        # Start FastAPI server
        from app.main import app
        uvicorn.run(app, host=host, port=port, reload=False, log_level="info")

    except KeyboardInterrupt:
        print("\n\n  Server terminated by user (Ctrl+C). Goodbye.")
    except Exception as e:
        import traceback
        print("\n" + "!" * 76)
        print("  AN UNEXPECTED ERROR OCCURRED DURING EXECUTION:")
        print("!" * 76)
        traceback.print_exc()
        print("\n  If the issue persists, please report to: drugdesignsynthesislab@gmail.com")
        input("\n  Press Enter to close this window...")

if __name__ == "__main__":
    main()
