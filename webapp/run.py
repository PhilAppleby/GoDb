from app import app
port=8085
print("Starting on port {0}".format(port))

if __name__ == "__main__":
  app.run(port=port)
