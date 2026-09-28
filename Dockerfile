FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV WINEPREFIX=/root/.wine
ENV WINEDEBUG=-all
ENV DISPLAY=:99
ENV PORT=7860

# Enable i386 multiarch and install wine32, Xvfb, python3, and mingw compiler
RUN dpkg --add-architecture i386 && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        wine32 \
        wine \
        xvfb \
        python3 \
        python3-pip \
        gcc-mingw-w64-i686 \
        git && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Pre-initialize Wine environment
RUN wineboot --init || true

# Install Python requirements
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy all project files including bundled models and CORALSEA.exe
COPY . .

# Build the native Win32 driver executable
RUN i686-w64-mingw32-gcc -O2 coral_runner.c -o coral_runner.exe

EXPOSE 7860

# Launch headless display server and FastAPI prediction engine
CMD ["sh", "-c", "Xvfb :99 -screen 0 1024x768x16 & python3 run_server.py"]
