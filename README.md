# Python project for calculations of **Close Earth Aproaches**
It was created as a university homework as a part of OZ0-A course.

## Currently hosted 
Should be up and running, see for your self:
[https://closeapproach-production.up.railway.app](https://closeapproach-production.up.railway.app) :point_up:

> If there is something wrong let me know!

## Setup virtual environment to **run localy**
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

### Uninstall
Simply remove the repository and uninstall the python@3.10 if you don't need it anymore.


## Issues
- [] Axes are not limited, it is issue when ploting hyperbolic orbit
- [] Appearance of the page is limited :laughing:
- [] Bodies added to the plot cannot be changed separately
