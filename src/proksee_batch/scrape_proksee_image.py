import glob
import os
import time

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager


def setup_browser(download_dir: str) -> webdriver.Chrome:
    """
    Sets up the Chrome browser for headless operation.
    """
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--window-size=1920,1080")

    os.makedirs(download_dir, exist_ok=True)
    chrome_options.add_experimental_option(
        "prefs",
        {
            "download.default_directory": download_dir,
            "download.prompt_for_download": False,
            "download.directory_upgrade": True,
            "safebrowsing.enabled": True,
        },
    )

    service = Service(ChromeDriverManager().install())
    return webdriver.Chrome(service=service, options=chrome_options)  # type: ignore


def download_data(browser: webdriver.Chrome, url: str) -> None:
    """
    Navigates to the given URL and performs the downloading actions.
    """
    browser.get(url)
    time.sleep(5)

    try:
        first_button = browser.find_element("id", "react-tabs-4")
        first_button.click()
        time.sleep(3)

        second_button = browser.find_element(
            "css selector", "img[title='Download SVG']"
        )
        second_button.click()
        time.sleep(5)
    finally:
        browser.quit()


def scrape_proksee_image(proksee_link_file: str, output_file: str) -> None:
    """
    Scrapes the Proksee image from the given Proksee link file and saves it to the given output file.
    """
    # Define a temporary directory
    download_dir = os.path.join(
        os.path.dirname(output_file),
        "temporary_" + proksee_link_file.rsplit("/", 1)[-1],
    )

    # Extract link from file.
    proksee_link = None
    with open(proksee_link_file) as file:
        proksee_link = file.read().strip()

    # Set up the browser and download the data.
    browser = setup_browser(download_dir)
    download_data(browser, proksee_link)

    # Check that the download was successful.
    if len(glob.glob(os.path.join(download_dir, "*.svg"))) == 0:
        raise Exception("No .svg file was downloaded.")

    # Move the downloaded .svg file from the temporary directory to the output path.
    svg_file = glob.glob(os.path.join(download_dir, "*.svg"))[0]
    os.rename(svg_file, output_file)

    # Remove the temporary directory.
    os.rmdir(download_dir)
