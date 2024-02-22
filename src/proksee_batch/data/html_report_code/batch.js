////////////////////////////////////////////////////////////////////////////////
// Initial setup
////////////////////////////////////////////////////////////////////////////////

// NOTE: tableData (global variable) is the main data object that contains all the genome data
// Adjust table data to include keys
const genomeData = tableData.genomes;
Object.keys(genomeData).forEach(key => {
    genomeData[key].key = key;
});

updateRunInfo();

// Current Sort Settings
let currentSort = { column: null, isDescending: false };
// This will hold a filtered version of the table data when a search is performed
let filteredData = {};
// Current search terms
let searchTerms = [];

// Create viewer
const cgv = new CGV.Viewer('#my-viewer', {height: 500});
window.cgv = cgv;

// Reset map and search
autoResizeMyViewer();
clearMap();
updateSearchCount();

// Initial call to generate table
generateGenomeList();

////////////////////////////////////////////////////////////////////////////////
// Main function to generate the genome list
////////////////////////////////////////////////////////////////////////////////

function generateGenomeList() {
    let data = (searchTerms.length > 0) ? filteredData : genomeData;
    data = sortData(data);

    console.log('Generating genome list:', data);
    const genomeList = document.getElementById('genome-list-container');
    genomeList.innerHTML = '';
    const items = Array.isArray(data) ? data : Object.entries(data);
    items.forEach(([key, item], index) => {
        const contigCount = item["Number of contigs"];
        const genomeName = highlightSearchTerms(item.Name, searchTerms);
        const genomeDescription = highlightSearchTerms(item.Description, searchTerms);
        const selected = (key == selectedMapKey()) ? 'selected' : '';
        const fileTags = getFileTags(item);
        genomeList.innerHTML += `
            <div class='genome-item ${selected}' data-key='${key}'>
                <div class='genome-left'>
                    <div class='genome-name'>${genomeName}</div>
                    <div class='genome-description'>${genomeDescription}</div>
                    <div class='genome-tags'>${fileTags}</div>
                </div>
                <div class='genome-right'>
                    <div class='genome-size'>${item["Total size"].toLocaleString()} bp</div>
                    <div class='genome-count'><span class='genome-count-number'>${contigCount}</span> ${(contigCount > 1) ? 'contigs' : 'contig'}</div>
                    <div class='genome-gc'>${(item["GC content"] * 100).toFixed(2)}% GC</div>
                </div>
            </div>
        `;
    });

    // Add event listeners to each genome item
    const genomeItems = document.querySelectorAll('.genome-item');
    genomeItems.forEach(item => {
        item.addEventListener('click', function(event) {
            // Could also use event.target.dataset.key
            const key = event.currentTarget.getAttribute('data-key');
            if (key) {
                console.log('Clicked on genome:', key);
                loadDataForKey(key);
            }
        });
    });
}

function getFileTags(genomeItem) {
    const tagMap = {
        'genbank': 'GenBank',
        'json': 'JSON',
        'blast': 'BLAST',
        'bed': 'BED',
        'vcf': 'VCF',
        'gff': 'GFF',
    };
    console.log("Item:", genomeItem)
    if (!genomeItem.Files) return '';
    const fileTags = Object.keys(genomeItem.Files).map(fileType => {
        const tagType = tagMap[fileType] || fileType;
        const fileCount = genomeItem.Files[fileType].length;
        const filesString = genomeItem.Files[fileType].join(', ');
        const tagSuffix = (fileCount > 1) ? `<span class='tag-suffix'> x${fileCount}</span>` : '';
        return (fileCount > 0) ? `<div class='tag tag-${fileType}' title='${filesString}'>${tagType}${tagSuffix}</div>` : '';
    });
    return fileTags.join('');
};

////////////////////////////////////////////////////////////////////////////////
// Sorting functions
////////////////////////////////////////////////////////////////////////////////

const sortButtons = document.querySelectorAll('.sort-btn');
const sortMap = {
    'sort-name': 'Name',
    'sort-length': 'Total size',
    'sort-contigs': 'Number of contigs',
    'sort-gc': 'GC content',
};
sortButtons.forEach(btn => {
    btn.addEventListener('click', function(event) {
        const target = event.target;
        const column = sortMap[target.id];
        updateSortSettings(column);
        adjustSortIndicators(target);
        generateGenomeList()
    });
});

// function getSortIndicator(column) {
//     return (currentSort.column === column) ? (currentSort.isDescending ? '↓' : '↑') : '';
// }

function adjustSortIndicators(target) {
    sortButtons.forEach(btn => {
        btn.classList.remove('sort-asc');
        btn.classList.remove('sort-desc');
    });
    if (target.id === 'sort-off') return;
    if (currentSort.isDescending) {
        target.classList.add('sort-desc');
    } else {
        target.classList.add('sort-asc');
    }
}

function updateSortSettings(column) {
    const modifier = currentSort.isDescending ? -1 : 1;

    if (currentSort.column !== column) currentSort.isDescending = false;
    else currentSort.isDescending = !currentSort.isDescending;

    currentSort.column = column;
}

function sortData(data) {
    if (!currentSort.column) return data;
    const column = currentSort.column;
    const modifier = currentSort.isDescending ? -1 : 1;

    const isNumeric = ['Total size', 'Number of contigs', 'GC content'].includes(column);
    let dataArray = Object.values(data);
    dataArray.sort((a, b) => {
        const valueA = isNumeric ? a[column] : a[column].toUpperCase();
        const valueB = isNumeric ? b[column] : b[column].toUpperCase();
        return (valueA < valueB ? -1 : valueA > valueB ? 1 : 0) * modifier;
    });

    // Convert the sorted array back to an object
    let sortedData = {};
    dataArray.forEach(item => {
        sortedData[item.key] = item;
    });

    return sortedData;
}


////////////////////////////////////////////////////////////////////////////////
// Code for generating Proksee project links
////////////////////////////////////////////////////////////////////////////////

function generateProkseeLink(element, sampleId) {
    loadScript(`${sampleId}.js`, function() {
        if (typeof window.json === 'undefined') {
            console.error('No data found for sample ID:', sampleId);
            alert('Failed to generate Proksee project: No data available for this sample.');
            return;
        }

        const data = { origin: 'proksee-batch', data: JSON.stringify(window.json) };
        const url = 'https://proksee.ca/api/v1/projects.json';
        fetch(url, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data)
        })
        .then(response => response.json())
        .then(data => {
            if (data?.status === 'success' && data?.url) {
                const link = document.createElement('a');
                link.textContent = 'Go to Proksee Project';
                link.href = data.url;
                link.className = 'generated-link';
                link.target = '_blank';
                element.parentNode.replaceChild(link, element);
            } else {
                console.error(`Failed to create Proksee project: ${data?.error}`);
            }
        })
        .catch(error => {
                console.error('Error:', error);
        });
    });
}

function loadScript(scriptUrl, callback) {
    const existingScript = document.querySelector(`script[src="${scriptUrl}"]`);
    if (existingScript) {
        existingScript.remove();
    }
    const script = document.createElement('script');
    script.src = scriptUrl;
    script.type = 'text/javascript';
    document.body.appendChild(script);
    script.onload = function() {
        console.log(`Script loaded: ${scriptUrl}`);
        if (typeof callback === "function") {
            callback();
        }
    };
    script.onerror = function() {
        console.error(`Error loading script: ${scriptUrl}`);
        alert(`Failed to load data for script: ${scriptUrl}`);
    };
}


////////////////////////////////////////////////////////////////////////////////
// Code for loading CGView map from relevant .js files when table rows are clicked.
////////////////////////////////////////////////////////////////////////////////

function loadDataForKey(key) {
    const myViewer = document.querySelector('#my-viewer');
    const mapName = document.querySelector('.map-genome-name');
    showMessage('<div class="spinner-container"><div class="spinner"></div>Loading...</div>');
    mapName.textContent = genomeData[key].Name;
    highlightSelectedGenome(key);
    var script = document.createElement('script');
    const dataPath = `data/${key}.js`;
    script.src = dataPath;
    script.onload = function() {
        console.log('Loading data file:', dataPath)
        setTimeout(() => {
            myViewer.dataset.key = key;
            cgv.io.loadJSON(window.json);
            // NOTE: resizing fixes an issue where the previous legend is not removed
            cgv.resize();
            cgv.draw();
            hideMessage();
        }, 1);
    };
    document.head.appendChild(script);
}

function selectedMapKey() {
    const myViewer = document.querySelector('#my-viewer');
    return myViewer.dataset.key;
}

function highlightSelectedGenome(key) {
    document.querySelectorAll('.genome-item').forEach(item => {
        item.classList.remove('selected');
    });
    document.querySelector(`.genome-item[data-key="${key}"]`).classList.add('selected');
}


////////////////////////////////////////////////////////////////////////////////
// Search Bar
////////////////////////////////////////////////////////////////////////////////

const searchInput = document.getElementById('search-input');
searchInput.addEventListener('input', function() {
    const searchCancel = document.getElementById('search-cancel');
    searchCancel.style.display = (searchInput.value !== '') ? 'block' : 'none';
    searchGenomes(searchInput.value);
    updateSearchCount()

});
const searchCancel = document.getElementById('search-cancel');
searchCancel.addEventListener('click', function() {
    searchInput.value = '';
    searchInput.focus();
    searchCancel.style.display = 'none';
    searchGenomes();
    updateSearchCount()
});

function searchGenomes(searchString = '') {
    searchTerms = searchString.trim().toLowerCase().split(/\s+/);
    const dataArray = Object.keys(genomeData).map(key => {
        const item = genomeData[key];
        item.key = key;
        return item;
    });
    console.log('Search terms:', searchTerms);
    if (searchTerms.length > 0) {
        const filteredArray = dataArray.filter(item => {
            itemSearchSpace = item.Name + ' ' + item.Description;
            return searchTerms.every(term => itemSearchSpace.toLowerCase().includes(term));
        });
        filteredData = {};
        filteredArray.forEach(item => {
            filteredData[item.key] = item;
        });
    } else {
        filteredData = genomeData;
    }
    // generateGenomeList(filteredData, searchTerms);
    generateGenomeList();
}

function updateSearchCount() {
    const searchInput = document.getElementById('search-input');
    const searchCount = document.getElementById('genome-count');
    totalGenomesCount = Object.keys(genomeData).length;
    filteredGenomesCount = Object.keys(filteredData).length;
    console.log('Total genomes:', totalGenomesCount);
    console.log('Search input:', searchInput.value);
    if (searchInput.value === '') {
        searchCount.textContent = `Genomes total: ${totalGenomesCount}`;
    } else {
        searchCount.textContent = `Genomes found: ${filteredGenomesCount} (of ${totalGenomesCount})`;
    }
}

function highlightSearchTerms(text, terms) {
    if (terms.length === 0) return text;
    const regex = new RegExp(`(${terms.join('|')})`, 'gi');
    return text.replace(regex, '<mark>$1</mark>');
}


////////////////////////////////////////////////////////////////////////////////
// Auto Resize and Drag Bar
////////////////////////////////////////////////////////////////////////////////

function mapResize() {
    const myViewer = document.querySelector('#my-viewer');
    const height = myViewer.clientHeight - 1;
    const width = myViewer.clientWidth - 1;
    cgv.resize(width, height);
}

function autoResizeMyViewer() {
  window.onresize = mapResize;
  window.onload = function () {
    setTimeout( () => {
      mapResize();
    }, 100);
  }
}

// Drag Bar
const resizeBar = document.querySelector('#resize-bar');
const listPanel = document.querySelector('.section-table');
const mapPanel = document.querySelector('.section-map');
resizeBar.addEventListener('mousedown', function(e) {
    e.preventDefault();
    document.addEventListener('mousemove', onMouseDragMove);
    document.addEventListener('mouseup', onMouseDragUp);
});

function onMouseDragMove(e) {
    const myViewer = document.querySelector('#my-viewer');
    const containerOffsetLeft = listPanel.offsetLeft;
    const newLeftWidth = e.clientX - containerOffsetLeft;
    listPanel.style.width = `${newLeftWidth}px`;
    mapPanel.style.width = `calc(100% - ${newLeftWidth}px)`;
    showMessage('Resizing...');
    myViewer.style.display = 'none';
}

function onMouseDragUp() {
    const myViewer = document.querySelector('#my-viewer');
    myViewer.style.display = 'block';
    // Only hide message if the map is visible
    const selectedKey = selectedMapKey();
    if (selectedKey) {
        hideMessage();
    } else {
        showMessage('Click on a genome to view map...');
    }
    mapResize();
    document.removeEventListener('mousemove', onMouseDragMove);
    document.removeEventListener('mouseup', onMouseDragUp);
}

////////////////////////////////////////////////////////////////////////////////
// Run Summary/Details
////////////////////////////////////////////////////////////////////////////////

function updateRunInfo() {
    const inputDirBase = tableData.input_dir.split('/').pop();
    const runName = tableData.run_name || inputDirBase;
    const runDate = tableData.run_date;
    const runVersion = tableData.version || '1.0.0missing';

    // Summary
    const runNameEl = document.querySelector('.run-name');
    const runDateEl = document.querySelector('.run-date');
    runNameEl.innerHTML = runName;
    runDateEl.innerHTML = runDate;
    // Details
    const detailVersionEl = document.querySelector('.batch-version');
    const detailDateEl = document.querySelector('.batch-run-date');
    const detailInputDirEl = document.querySelector('.batch-input-dir');
    detailVersionEl.innerHTML = runVersion;
    detailDateEl.innerHTML = runDate;
    detailInputDirEl.innerHTML = tableData.input_dir;

    const runSummaryEl = document.querySelector('.run-summary');
    const runDetailsEl = document.querySelector('.run-details');
    runSummaryEl.addEventListener('click', function() {
        if (runSummaryEl.classList.contains('active')) {
            runDetailsEl.classList.remove('active');
            runSummaryEl.classList.remove('active');
        } else {
            runDetailsEl.classList.add('active');
            runSummaryEl.classList.add('active');
        }
    });
}


////////////////////////////////////////////////////////////////////////////////
// Helpers: show/hide messages, clear map
////////////////////////////////////////////////////////////////////////////////

function clearMap() {
    cgv.io.loadJSON({cgview: {version: "1.6.0"}}); // Clear the map
    cgv.draw();
    showMessage('Click on a genome to view map...');
}

function showMessage(message) {
    const messageContainer = document.querySelector('#my-viewer-message');
    messageContainer.style.opacity = 1;
    messageContainer.style.display = 'flex';
    messageContainer.innerHTML = `<div id='message'>${message}</div>`;
}
function hideMessage() {
    const messageContainer = document.querySelector('#my-viewer-message');
    messageContainer.style.opacity = 0;
    setTimeout(() => {
        messageContainer.style.display = 'none';
    }, 1000); // This should match the transition duration in the CSS (or be a little bigger)
}
