# virus-shiny-app

## Docker

### Build docker image
```
docker build -t <image-name, e.g. ljelonek/virus:latest> .
```

### Run docker image
```
docker run -it --rm -p <localport>:3838 <image>
```

Explanation:
* `-it` bind container input and output to current terminal. Use `-d` if you
  want to run it a a daemon
* `--rm` delete container when it stops
* `-p <host-port>:<container-port>` Allows access to the container-app via
  `http://localhost:<host-port>`. The shiny app is configured to run with port `3838`
  in the container

Other useful options:
* `-v <host-path>:<container-path>` mounts a directory from the host computer
  into the container
