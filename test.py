#!/usr/bin/env python3
"""
test_upload.py

A simple test script to verify GitHub upload and Python execution.
"""

import sys
from datetime import datetime


def main():
    print("✅ Test script executed successfully!")
    print("-" * 50)

    print("Python version:")
    print(sys.version)

    print("\nCurrent date and time:")
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    print("\nThis script was uploaded to GitHub for testing purposes.")
    print("If you can see this output, your upload worked!")


if __name__ == "__main__":
    main()
