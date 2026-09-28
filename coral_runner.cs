using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Diagnostics;
using System.Runtime.InteropServices;

class CoralRunner
{
    private const uint BM_CLICK = 0x00F5;

    [DllImport("user32.dll", SetLastError = true)]
    static extern bool EnumWindows(EnumWindowsProc lpEnumFunc, IntPtr lParam);

    [DllImport("user32.dll", SetLastError = true)]
    static extern bool EnumChildWindows(IntPtr hWndParent, EnumChildProc lpEnumFunc, IntPtr lParam);

    [DllImport("user32.dll", SetLastError = true)]
    static extern uint GetWindowThreadProcessId(IntPtr hWnd, out uint lpdwProcessId);

    [DllImport("user32.dll", CharSet = CharSet.Auto, SetLastError = true)]
    static extern int GetClassName(IntPtr hWnd, StringBuilder lpClassName, int nMaxCount);

    [DllImport("user32.dll", CharSet = CharSet.Auto, SetLastError = true)]
    static extern int GetWindowText(IntPtr hWnd, StringBuilder lpString, int nMaxCount);

    [DllImport("user32.dll", CharSet = CharSet.Auto)]
    static extern IntPtr PostMessage(IntPtr hWnd, uint Msg, IntPtr wParam, IntPtr lParam);

    delegate bool EnumWindowsProc(IntPtr hWnd, IntPtr lParam);
    delegate bool EnumChildProc(IntPtr hWnd, IntPtr lParam);

    static int Main(string[] args)
    {
        string runDir = (args.Length >= 1 && args[0] != ".") ? Path.GetFullPath(args[0]) : Directory.GetCurrentDirectory();
        int timeoutSec = args.Length >= 3 ? int.Parse(args[2]) : 25;

        Directory.SetCurrentDirectory(runDir);

        string exePath = "CORALSEA.exe";
        string outPath = "ListModel.txt";

        if (!File.Exists(exePath))
        {
            Console.Error.WriteLine("Executable not found: " + exePath);
            return 2;
        }

        if (File.Exists(outPath))
        {
            try { File.Delete(outPath); } catch { }
        }

        ProcessStartInfo psi = new ProcessStartInfo(exePath);
        psi.WorkingDirectory = runDir;
        psi.UseShellExecute = false;

        Process p = Process.Start(psi);
        if (p == null)
        {
            Console.Error.WriteLine("Failed to start CORALSEA.exe");
            return 3;
        }

        try
        {
            IntPtr targetHwnd = IntPtr.Zero;
            Stopwatch sw = Stopwatch.StartNew();
            while (sw.ElapsedMilliseconds < 8000 && targetHwnd == IntPtr.Zero)
            {
                EnumWindows((hWnd, lp) =>
                {
                    uint pid;
                    GetWindowThreadProcessId(hWnd, out pid);
                    if (pid == p.Id)
                    {
                        StringBuilder sb = new StringBuilder(256);
                        GetClassName(hWnd, sb, 256);
                        if (sb.ToString() == "TForm1")
                        {
                            targetHwnd = hWnd;
                            return false;
                        }
                    }
                    return true;
                }, IntPtr.Zero);

                if (targetHwnd == IntPtr.Zero) Thread.Sleep(150);
            }

            if (targetHwnd == IntPtr.Zero)
            {
                Console.Error.WriteLine("Timed out waiting for TForm1");
                return 4;
            }

            IntPtr btnLoad = IntPtr.Zero;
            sw.Restart();
            while (sw.ElapsedMilliseconds < 4000 && btnLoad == IntPtr.Zero)
            {
                btnLoad = FindChild(targetHwnd, "Load method");
                if (btnLoad == IntPtr.Zero) Thread.Sleep(150);
            }

            if (btnLoad != IntPtr.Zero)
            {
                PostMessage(btnLoad, BM_CLICK, IntPtr.Zero, IntPtr.Zero);
                Thread.Sleep(350);
            }

            IntPtr btnImport = FindChild(targetHwnd, "Import of current model");
            if (btnImport != IntPtr.Zero)
            {
                PostMessage(btnImport, BM_CLICK, IntPtr.Zero, IntPtr.Zero);
                Thread.Sleep(400);
            }

            IntPtr btnCalc = FindChild(targetHwnd, "Calculation model for a list of SMILES");
            if (btnCalc != IntPtr.Zero)
            {
                PostMessage(btnCalc, BM_CLICK, IntPtr.Zero, IntPtr.Zero);
            }

            sw.Restart();
            bool success = false;
            while (sw.ElapsedMilliseconds < (timeoutSec * 1000))
            {
                Thread.Sleep(150);
                if (File.Exists(outPath) && new FileInfo(outPath).Length > 0)
                {
                    success = true;
                    Thread.Sleep(150);
                    break;
                }
            }

            if (success)
            {
                Console.WriteLine("SUCCESS: ListModel.txt generated.");
                return 0;
            }
            else
            {
                Console.Error.WriteLine("Timed out waiting for ListModel.txt");
                return 5;
            }
        }
        finally
        {
            try
            {
                if (!p.HasExited)
                {
                    p.Kill();
                    p.WaitForExit(1000);
                }
            }
            catch { }
        }
    }

    static IntPtr FindChild(IntPtr parent, string prefix)
    {
        IntPtr result = IntPtr.Zero;
        EnumChildWindows(parent, (hWnd, lp) =>
        {
            StringBuilder sb = new StringBuilder(256);
            GetWindowText(hWnd, sb, 256);
            if (sb.ToString().StartsWith(prefix))
            {
                result = hWnd;
                return false;
            }
            return true;
        }, IntPtr.Zero);
        return result;
    }
}
