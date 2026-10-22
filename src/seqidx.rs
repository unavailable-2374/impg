use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SequenceIndex {
    pub name_to_id: FxHashMap<String, u32>,
    id_to_name: FxHashMap<u32, String>,
    id_to_len: FxHashMap<u32, usize>,
    next_id: u32,
}

impl SequenceIndex {
    /// Cache accounting estimate, including both owned copies of each name.
    pub(crate) fn estimated_heap_bytes(&self) -> usize {
        self.name_to_id.capacity() * (std::mem::size_of::<(String, u32)>() + 2)
            + self.id_to_name.capacity() * (std::mem::size_of::<(u32, String)>() + 2)
            + self.id_to_len.capacity() * (std::mem::size_of::<(u32, usize)>() + 2)
            + self.name_to_id.keys().map(|s| s.capacity()).sum::<usize>()
            + self.id_to_name.values().map(|s| s.capacity()).sum::<usize>()
    }
    pub fn new() -> Self {
        SequenceIndex {
            name_to_id: FxHashMap::default(),
            id_to_name: FxHashMap::default(),
            id_to_len: FxHashMap::default(),
            next_id: 0,
        }
    }

    pub fn get_or_insert_id(&mut self, name: &str, length: Option<usize>) -> u32 {
        let id = *self.name_to_id.entry(name.to_owned()).or_insert_with(|| {
            let id = self.next_id;
            self.id_to_name.insert(id, name.to_owned());
            self.next_id += 1;
            id
        });

        if let Some(len) = length {
            self.id_to_len.entry(id).or_insert(len);
        }

        id
    }

    pub fn get_id(&self, name: &str) -> Option<u32> {
        self.name_to_id.get(name).copied()
    }

    pub fn get_name(&self, id: u32) -> Option<&str> {
        self.id_to_name.get(&id).map(|s| s.as_str())
    }

    pub fn get_len_from_id(&self, id: u32) -> Option<usize> {
        self.id_to_len.get(&id).copied()
    }

    pub fn is_empty(&self) -> bool {
        self.name_to_id.is_empty()
    }

    pub fn len(&self) -> usize {
        self.name_to_id.len()
    }
}
