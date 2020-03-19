if [ "$(docker ps -aq -f name=nanome-hydrogens)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f nanome-hydrogens
fi

docker run -d \
--name nanome-hydrogens \
--restart unless-stopped \
-e ARGS="$*" \
nanome-hydrogens
