import os
import sys
import webbrowser
import uvicorn

if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

def main():
    port = int(os.environ.get("PORT", 8000))
    host = os.environ.get("HOST", "0.0.0.0")
    display_url = f"http://localhost:{port}"

    print("=" * 76)
    print("  Multi-Kinase 2D-QSAR Prediction Platform (CORALSEA Monte Carlo)")
    print("  Drug Design Synthesis Lab")
    print("  Department of Pharmaceutical Sciences and Drug Research")
    print("  Punjabi University, Patiala, Punjab, India")
    print("  Contact: drugdesignsynthesislab@gmail.com")
    print("=" * 76)
    print(f"\n  Starting Web Portal at: {display_url}")
    print("  Active Kinase Engine: NF-κB Inducing Kinase (NIK / MAP3K14)")
    print("  Press Ctrl+C to terminate the server.\n")

    # Launch browser automatically
    try:
        webbrowser.open(display_url)
    except Exception:
        pass

    # Start FastAPI server
    uvicorn.run("app.main:app", host=host, port=port, reload=False, log_level="info")

if __name__ == "__main__":
    main()
