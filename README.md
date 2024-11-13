## Setup virtual environment
To run this project it is recommended to create virtual environment and install required modules from *requirements.txt*. This procedure will ensures the separation of your system Python and this project. After installation the virtual will contain its own instance of specific Python version and required packages.

### Python3 installation
On Windows, download and add to yout PATH variable:
[https://www.python.org/ftp/python/3.10.0/python-3.10.0-amd64.exe](https://www.python.org/ftp/python/3.10.0/python-3.10.0-amd64.exe)

On UNIX system use prefered package manager, example:
```shell
    brew install python@3.10
```
### Python dependencies installation
In the project directory run
```shell
    python3 -m venv venv
```
On Windows:
```shell
    source venv/bin/activate
```
On UNIX:
```shell
    venv\Scripts\activate
```
You are within virtual environment, now install dependencies
```shell
    pip3 install -r requirements.txt
```

### Run code 
```shell
    python3 ./app.py
```

### If nothing happens
The application is running on [http://127.0.0.1:5000](http://127.0.0.1:5000), try to look at it on your browser.
