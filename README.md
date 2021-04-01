# Nanome - Hydrogens

A Nanome Integration Plugin to add/remove hydrogens to/from structures.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

## Usage

To run Hydrogens in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address> [optional args]
```

## Development

To run Hydrogens with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> [optional args]
```

## License

MIT
